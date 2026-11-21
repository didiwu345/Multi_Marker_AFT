# Multi_Marker_AFT

## Introduction to codebase
`Multi_Marker_AFT.r` is the main function that implements multi-marker association and interaction test based on the AFT model.

To successfully implement the test, you will need following R package dependencies:
* MASS
* CompQuadForm
* survival
* matrixStats
* matrixcalc
* RcppArmadillo
* Rcpp
* Matrix
* lbaft
* coxKM

`sample_data` folder contains essential sample datasets to conduct an example association test, without considering genetic heterogeneity, with covariate adjustment and left truncation.

`utils` folder contains utility C++ functions `aft_comprisk.cpp` and `aft_comprisk_smalln.cpp` that will be used in main function `Multi_Marker_AFT.r`.

`Missing_genotype_imputation` folder has a sample script to perform missing genotype imputation with IMPUTE2 software. The best practices using IMPUTE2 software can be found in:
[Getting started with IMPUTE2](https://mathgen.stats.ox.ac.uk/impute/impute_v2.html#getting_started). A sample code that executes IMPUTE2 software (v2.3.2) code is provided in `Missing_genotype_imputation/imputation.txt`.

## Multi_Marker_AFT function input and output
###Input:
* `G`:   Genetic variant matrix
* `H`:   Additional genetic variant matrix (G-G) or environment matrix (G-E) (Default: NULL; Only needed in G-G/G-E interaction). When `H` is given, an interaction test is performed: centered `G` and `H` main effects are added to `adj_cov` in the null model, and only the G-H interaction effect is tested.
* `GxH`: Interaction terms as input matrix (Default: NULL). If NULL and `H` is given, `GxH` is built from all pairwise products of the columns of centered `G` and centered `H`.
* `adj_cov`:   Adjustment covariate matrix (Default: NULL; required for main-effect test)
* `het_cov`:   Covariate matrix indicating genetic heterogeneity
* `kernel_G`: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') or Kernel matrix. Main-effect test: kernel on `G`. Interaction test with `kernel_H=NULL`: kernel on `GxH` ('cp' or 'gaussian' recommended; 'ibs' assumes 0/1/2 genotypes). Interaction test with `kernel_H` given: kernel on `G`, and the interaction kernel is `K_G * K_H` (element-wise).
* `kernel_H`: Kernel type   (included: 'ibs', 'cp', 'gaussian', 'identity') for H or Kernel matrix (Default: NULL; Only used in G-G/G-E interaction). With `kernel_G='cp'` and `kernel_H='cp'`, `K_G * K_H` equals the 'cp' kernel on `GxH`.
* `kernel_het`: Kernel type (included: 'ibs', 'cp', 'gaussian', 'identity') for het_cov or Kernel matrix (Default: NULL; Only needed under genetic heterogeneity)
* `smalln_ind`: 'smalln_adj' or 'no_smalln_adj' (Default)
* `trunct`: truncation time, default is NULL --> no truncation time (OR all truncation time is 0)
* `survt`: event or censor time
* `status`: Delta (0 or 1)
* `BB`: Number of iterations for W, S matrix approximation (Default: 500); Increase BB increases approximation accuracy, but also increases computing time.

###Output:
* `pval`: p-value


## Example of Multi_Marker_AFT method implementation

### Step1: Clone the repo
```bash
cd ~
git clone https://github.com/didiwu345/Multi_Marker_AFT.git
cd Multi_Marker_AFT
```

### Step 2: Install required R packages
Make sure you follow the [instruction](https://cran.r-project.org/bin/macosx/) to install R-4.3.3, as well as downloading and installing XQuartz and GNU Fortran. Then install the required CRAN or non-CRAN R packages. For non-CRAN packages, the corresponding source package is provided in the repo.
```
install.packages('quadprog')
install.packages('broom')
install.packages('magrittr')
install.packages('./BB_2019.10-1.tar.gz', repos = NULL, type="source")
install.packages('./geepack_1.3.10.tar.gz', repos = NULL, type="source")
install.packages('./lbaft_0.0-1.tar.gz', repos = NULL, type="source")
install.packages('./coxKM_0.3.tar.gz', repos = NULL, type="source")
```

### Step3: Sample implementation
This code example implements the `Multi_Marker_AFT` method in detecting genetic association effects, with covariate adjustment and left truncation, without considering confounding effect and genetic heterogeneity. 
The sample genotype dataset has sample size `N=500`, and SNP-set size `p=10`. Since `N >> p`, small sample correction is turned off by setting `smalln_ind='no_smalln_adj'`. `IBS` kernel is used to measure genetic similarity by setting `kernel_G = 'ibs'`.

```bash
suppressMessages(library(MASS))
suppressMessages(library(CompQuadForm))
suppressMessages(library(survival))
suppressMessages(library(matrixStats))
suppressMessages(library(matrixcalc))
suppressMessages(library(RcppArmadillo))
suppressMessages(library(Rcpp))
suppressMessages(library(Matrix))
suppressMessages(library(lbaft))
suppressMessages(library(coxKM))

sourceCpp('utils/aft_comprisk_smalln.cpp')
ibskm = function(X){ return(1 - as.matrix(dist(X, method = "manhattan"))/(2*ncol(X))) }  ### IBS kenrel: Input genotype matrix
cpkm = function(X){ return(tcrossprod(X)) }                                              ### Cross-product kernel: Input genotype matrix
simkm = function(X){ return(exp(-as.matrix(dist(X,method='euclidean')^2)/ncol(X))) }     ### Gaussian kernel: Input covariate matrix
idtkm = function(X){ return(0 + outer(as.numeric(X), as.numeric(X), "==")) }             ### For 1-column vector
plkm = function(X){ return(ncol(X) + tcrossprod(X)) }                                    ### Product linear kernel

sample_data = readRDS('sample_data/sample_data.rds')

G = sample_data[['genotype']]
X = sample_data[['adj_cov']]
status = sample_data[['status']]
survt = sample_data[['survt']]
trunct = sample_data[['trunct']]

source('Multi_Marker_AFT.r')

pval = Multi_Marker_AFT(G=G, H=NULL, GxH=NULL, het_cov=NULL, adj_cov=X, kernel_G='ibs', kernel_H=NULL, kernel_het=NULL, smalln_ind='smalln_adj', trunct=trunct, survt=survt, status=status, BB=500)
```

### Step4: G-G / G-E interaction test
Pass the second SNP-set (G-G) or environment matrix (G-E) as `H`. The G and H main effects are adjusted in the null model automatically, so `adj_cov` only needs the other covariates.
```bash
#-- G-G interaction: cross-product kernel on GxH (built internally from centered G and H)
G1 = G[,1:5]
G2 = G[,6:10]
pval_GG = Multi_Marker_AFT(G=G1, H=G2, adj_cov=X, kernel_G='cp', smalln_ind='no_smalln_adj', trunct=trunct, survt=survt, status=status, BB=500)

#-- G-E interaction: IBS kernel on G times identity kernel on a binary environment E
E = as.matrix(X[,1])
pval_GE = Multi_Marker_AFT(G=G, H=E, adj_cov=X[,-1,drop=FALSE], kernel_G='ibs', kernel_H='identity', smalln_ind='no_smalln_adj', trunct=trunct, survt=survt, status=status, BB=500)
```

