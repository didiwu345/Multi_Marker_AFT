library(MASS)
library(CompQuadForm)
library(survival)
library(matrixStats)
library(matrixcalc)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(lbaft)
library(coxKM)

ibskm = function(X){ return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) }  ### IBS kenrel: Input genotype matrix
cpkm = function(X){ return(tcrossprod(X)) }                                              ### Cross-product kernel: Input genotype matrix
simkm = function(X){ return(exp(-as.matrix(dist(X,method='euclidean')^2)/ncol(X))) }     ### Gaussian kernel: Input covariate matrix
idtkm = function(X){ return(0 + outer(as.numeric(X), as.numeric(X), "==")) }             ### Identity kernel for 1-column vector
plkm = function(X){ return(ncol(X) + tcrossprod(X)) }                                    ### Product linear kernel

### Build kernel matrix from a kernel type or return a user-supplied kernel matrix
get_kernel = function(kernel, X, arg_name){
  if(is.matrix(kernel)){
    return(kernel)
  } else if(is.null(X)){
    stop(paste0(arg_name, " is given as a kernel type but its input matrix is NULL"))
  } else if(kernel=='ibs'){
    return(ibskm(X))
  } else if(kernel=='cp'){
    return(cpkm(X))
  } else if(kernel=='gaussian'){
    return(simkm(X))
  } else if(kernel=='identity'){
    if(ncol(as.matrix(X)) != 1){ stop(paste0("'identity' kernel for ", arg_name, " requires a 1-column input")) }
    return(idtkm(X))
  } else {
    stop(paste0(arg_name, " should be either a matrix, or one of 'ibs','cp','gaussian','identity'"))
  }
}

### Load C++ utilities; both files define the same function names, so re-source only when the mode changes
load_cpp = function(cpp_file){
  if(!identical(getOption('mmaft_cpp_loaded'), cpp_file)){
    sourceCpp(cpp_file)
    options(mmaft_cpp_loaded = cpp_file)
  }
}

### Center columns of a matrix
center_cols = function(X){
  X = as.matrix(X)
  return(sweep(X, 2, colMeans(X)))
}

### Obtain p-value
Multi_Marker_AFT = function(G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=NULL, kernel_G, kernel_H=NULL, kernel_het=NULL, smalln_ind='no_smalln_adj',trunct=NULL,survt,status,BB=500){
  # Summary: Function implement HWV_AFT association and interaction test on provided inputs and outputs p-value
  # Input:
  #     G:   Genetic variant matrix
  #     H:   Additional genetic variant matrix (G-G) or environment matrix (G-E) (Default: NULL; Only needed in G-G/G-E interaction)
  #          When H is given, an interaction test is performed: G and H main effects are added to adj_cov in the null model,
  #          and only the G-H interaction effect is tested.
  #     GxH: Interaction terms as input matrix (Default: NULL). If NULL and H is given, GxH is built from all pairwise
  #          products of the columns of centered G and centered H.
  #     adj_cov:   Adjustment covariate matrix (Default: NULL; required for main-effect test)
  #     het_cov:   Covariate matrix indicating genetic heterogeneity
  #     kernel_G: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') or Kernel matrix.
  #               Main-effect test: kernel on G.
  #               Interaction test: kernel on G, and the interaction kernel is K_G * K_H (element-wise).
  #               If GxH is supplied and kernel_H=NULL: kernel on GxH instead ('cp' or 'gaussian' recommended).
  #               If kernel_G is a matrix and kernel_H=NULL: kernel_G is used as the full interaction kernel.
  #     kernel_H: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') for H or Kernel matrix
  #               (Default: NULL --> same kernel type as kernel_G, unless GxH is supplied; Only used in G-G/G-E interaction)
  #               With kernel_G='cp' and kernel_H='cp', K_G * K_H equals the 'cp' kernel on GxH built from centered G and H.
  #     kernel_het: Kernel type (included: 'ibs', 'cp', 'gaussian', 'identity') for het_cov or Kernel matrix (Default: NULL; Only needed under genetic heterogeneity)
  #     smalln_ind: 'smalln_adj' or 'no_smalln_adj' (Default)
  #     trunct: truncation time, default is NULL --> no truncation time (OR all truncation time is 0)
  #     survt: event or censor time
  #     status: Delta (0 or 1)
  #     BB: Number of iterations for W, S matrix approximation (Default: 500)
  # Output:
  #     pval: p-value

  G = as.matrix(G)
  n = nrow(G)
  if(is.null(trunct)){ trunct = rep(0, n) }

  #-- Interaction (G-G / G-E) setup
  interaction = !is.null(H) || !is.null(GxH)
  if(interaction){
    if(is.null(H)){
      stop("H must be provided for G-G/G-E interaction test, so that G and H main effects are adjusted in the null model")
    }
    Gc = center_cols(G)
    Hc = center_cols(H)
    #-- Without user-supplied GxH, the same kernel type is applied to G and H: K = K_G(G) * K_G(H)
    if(is.null(kernel_H) && is.null(GxH) && !is.matrix(kernel_G)){ kernel_H = kernel_G }
    if(is.null(GxH)){
      GxH = NULL
      for(i in 1:ncol(Hc)){ GxH = cbind(GxH, Gc*Hc[,i]) }
    } else {
      GxH = as.matrix(GxH)
    }
    #-- G and H main effects go into the null model
    adj_cov = cbind(Gc, Hc, adj_cov)
    #-- Drop covariates duplicating a main effect (e.g. E passed in both H and adj_cov), which makes the null model singular
    adj_cor = suppressWarnings(abs(cor(adj_cov)))
    adj_cor[is.na(adj_cor) | lower.tri(adj_cor, diag=TRUE)] = 0
    dup_cols = colSums(adj_cor > 1 - 1e-8) > 0
    if(any(dup_cols)){
      warning(paste(sum(dup_cols), "adj_cov column(s) duplicate a G/H main effect and were dropped from the null model"))
      adj_cov = adj_cov[, !dup_cols, drop=FALSE]
    }
  }

  if(is.null(adj_cov)){
    stop("adj_cov is required: the null AFT model needs at least one adjustment covariate")
  }
  adj_cov = as.matrix(adj_cov)

  #-- Similarity of tested effect
  if(!interaction){
    Gsim = get_kernel(kernel_G, G, 'kernel_G')
  } else if(is.null(kernel_H)){
    if(identical(kernel_G, 'ibs')){
      warning("'ibs' kernel assumes 0/1/2 genotypes and is not appropriate for interaction terms GxH; consider 'cp' or 'gaussian'")
    }
    Gsim = get_kernel(kernel_G, GxH, 'kernel_G')
  } else {
    Gsim = get_kernel(kernel_G, Gc, 'kernel_G') * get_kernel(kernel_H, Hc, 'kernel_H')
  }

  #-- Define K matrix
  if(is.null(kernel_het)){
    K = Gsim
  } else {
    Fsim = get_kernel(kernel_het, het_cov, 'kernel_het')
    K = Gsim * (1 + Fsim)
  }

  Krank = tryCatch(rankMatrix(K)[1],
                   error=function(e){ stop(paste("Rank computation of kernel matrix K failed:", conditionMessage(e))) },
                   warning=function(e){ stop(paste("Rank computation of kernel matrix K failed:", conditionMessage(e))) })

  K_eigen = eigen(K)
  K_evec = as.matrix(K_eigen$vectors)[,1:Krank]
  K_eval = abs(K_eigen$values[1:Krank])
  if(Krank==1){
    E1 = matrix(sqrt(K_eval)*K_evec, nrow=1)
  } else {
    E1 = diag(sqrt(K_eval)) %*% t(K_evec)
  }

  fit = suppressWarnings(fit.GP(Surv(as.matrix(trunct), as.matrix(survt), as.matrix(status)) ~ adj_cov - 1, vMethod='M1'))
  beta_est = as.matrix(fit$beta)

  #-- All time points (ordered)
  tall = as.numeric(sort(unlist(log(survt) - adj_cov %*% beta_est), decreasing = F))

  if(smalln_ind=='no_smalln_adj'){

    load_cpp('utils/aft_comprisk.cpp')

    W_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    Q_all <- Qall_cpp(E1, W_mats, adj_cov, trunct, survt, tall, status, beta_est)
    B_est = matrix(0, nrow=Krank, ncol=ncol(adj_cov))
    for(i in 1:Krank){
      fitlm = lm(Q_all[,i]/sqrt(n) ~ 0 + W_mats)
      B_est[i,] = as.numeric(fitlm$coefficients)
    }

    S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    U_all <- Uall_cpp(S_mats, adj_cov, trunct, survt, tall, status, beta_est)
    A_est = matrix(0,nrow=ncol(adj_cov),ncol=ncol(adj_cov))
    for(i in 1:ncol(adj_cov)){
      fitlm = lm(U_all[,i]*sqrt(n) ~ 0 + S_mats)
      A_est[i,] = as.numeric(fitlm$coefficients)
    }

    #-- Covariance matrix of E1M
    E1M_cov <- R1M_cov_cpp(E1, B_est, A_est, adj_cov, trunct, survt, tall, status, beta_est)
    #-- Eigen decomposition
    Chi_coef = abs(eigen(E1M_cov, symmetric=TRUE, only.values=TRUE)$values)
    #-- Test stat and p-value
    E1M = E1 %*% M_tild_cpp(adj_cov, trunct, survt, tall, status, beta_est)
    test_stat = c(crossprod(E1M))
    pval = davies(test_stat, Chi_coef)$Qq

  } else if(smalln_ind=='smalln_adj'){

    load_cpp('utils/aft_comprisk_smalln.cpp')

    cause1event = status
    M_tild = M_tild_cpp(adj_cov,trunct,survt,tall,cause1event,beta_est)
    q_star = c((t(M_tild) %*% K %*% M_tild) / crossprod(M_tild))
    Knew = K - q_star*diag(n)
    Knew_rank = rankMatrix(Knew)[1]
    Knew_eigen = eigen(Knew)
    R1 = as.matrix(Knew_eigen$vectors[,1:Knew_rank])
    Knew_eval = Knew_eigen$values[1:Knew_rank]

    W_mats <- mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    Q_all <- Qall_cpp(t(R1),W_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)
    B_est = matrix(0,nrow=Knew_rank,ncol=ncol(adj_cov))
    for(i in 1:Knew_rank){
      fitlm = lm(Q_all[,i]/sqrt(n) ~ 0+W_mats)
      B_m = as.numeric(fitlm$coefficients)
      B_est[i,] = B_m
    }

    S_mats = mvrnorm(BB, mu=rep(0,ncol(adj_cov)), Sigma=diag(ncol(adj_cov)))
    U_all = Uall_cpp(S_mats,adj_cov,trunct,survt,tall,cause1event,beta_est)
    A_est = matrix(0,nrow=ncol(adj_cov),ncol=ncol(adj_cov))
    for(i in 1:ncol(adj_cov)){
      fitlm = lm(U_all[,i]*sqrt(n) ~ 0+S_mats)
      A_p = as.numeric(fitlm$coefficients)
      A_est[i,] = A_p
    }

    #-- Covariance matrix of R1M
    R1M_cov = R1M_cov_cpp(t(R1),B_est,A_est,adj_cov,trunct,survt,tall,cause1event,beta_est)

    #-- R1M covariance Eigen decomposition
    R1M_cov_eig = eigen(R1M_cov,symmetric=TRUE)
    R1M_cov_evec = as.matrix(R1M_cov_eig$vectors)
    R1M_cov_eval = abs(R1M_cov_eig$values)

    #-- Get Chi-square distribution coefficients
    Amat = diag(sqrt(R1M_cov_eval)) %*% t(R1M_cov_evec) %*% diag(Knew_eval) %*% R1M_cov_evec %*% diag(sqrt(R1M_cov_eval))
    Chi_coef = eigen(Amat,symmetric=T,only.values=T)$values

    #-- Test stat and p-value
    pval = davies(0, Chi_coef)$Qq

  } else {
    stop("smalln_ind should be either 'smalln_adj' or 'no_smalln_adj'")
  }
  return(pval)
}

