
# lfmm
The R package **lfmm** provides function to estimate unobserved confounding
factors in large association studies. Considering that we observed n samples
of p explained variables, put in a matrix $Y$, and n sample of a variable of
interest. We want to find out which explained variables are associated with the
variable of interest. The variables of interest and other useful covariates for
the model are put in a matrix X. Let considering the factor augmented linear
model 

$$
Y = U V^T + X B^T + E,
$$

where

- Y is an n by p matrix of explained with n samples and p variables,
- X is an n by d matrix of the variable of interest for the association with
  other covariates useful for the model, 
- B is a d by p matrix of unobserved coefficients for the observed
  covariates
- U is an n by K matrix of latent factors,
- V is a K by p matrix of coefficients for the latent factord, 
- and E is an n by p matrix of residual error.

The R package **lfmm** provides function to estimate the U and V matrices
(`lfmm_ridge` and `lfmm_lasso`). After estimating confounding, we propose a
function to perform the association test by accounting for latent factors
(`lfmm_test`).

## Installation

Install the latest version from github (requires [devtools](https://github.com/hadley/devtools)):
```R
# install.packages("devtools")
devtools::install_github("bcm-uga/lfmm")
```

