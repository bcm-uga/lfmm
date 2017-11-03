
## latent factor mixed models: lfmm
The R package **lfmm** implements latent factor mixed models providing confounder 
adjusment for large-scale association studies. The model considers n samples of p response 
variables, stored in an n x p matrix, and one or several variables of interest (e.g., 
phenotypes, ecological variables). The package enables testing which response variables 
are associated with the variables of interest while removing confounding variation due to sample stratification or batch effects. The lfmm algorithms are fast and efficient and are supported by statistical
theory (oracle efficiency, identifiability, optimality of solutions).

The R package **lfmm** provides two functions to estimate K latent confounders (or factors):
`lfmm_ridge` and `lfmm_lasso`. Those functions are based on algorithms solving regularized
least squares problems. After estimation of the confounding factors, association testing
is performed by using the function `lfmm_test` and chi-square tests.

## Installation

Installing the latest version from github requires [devtools](https://github.com/hadley/devtools):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/lfmm")
```

