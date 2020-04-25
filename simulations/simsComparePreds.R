index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_compare_preds"
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, paste0("run-", index, ".rds"))
if(file.exists(output_file)){
    quit('no')
}

library(gtools)
library(here)
library(tidyverse)
library(compositions)
#source(here("scripts", "compRegression.R"))
source(file.path("..", "scripts", "compRegression.R"))
source(file.path("..", "scripts", "multFracReg.R"))


### KlD function
kld <- function(true, est) {
    d <- 0
    for(i in 1:length(true)) {
        if(true[i] == 0) {
            d <- d 
        } else {
            d <- d + true[i] * log(true[i] / est[i])
        }
    }
    return(d)
}

### Average KLD for compositional rows
kldMat <- function(trueMat, estMat) {
    d <- 0
    N <- nrow(trueMat)
    for(i in 1:N) {
        d <- d + kld(trueMat[i,], estMat[i,])
    }
    return(d / N)
}

C <- 3
### Create list of M for direct regression
### First choice of M is from education data set
M1 <- matrix(c(.91, .05, .04,
               0, .91, .09,
               0, .14, .86),
             nrow = C, byrow = TRUE)
### Second choice is from microscope data set
M2 <- matrix(c(.97, .03, 0,
               0, 1, 0,
              0, .04, .96),
             nrow = C, byrow = TRUE)

M_list <- list(M1, M2)

### List of betas for fractional logit model
beta_logit_list <- list(
    matrix(c(.06, .19, 0,
             -.78, .41,0,
             -.71, -.71, 0),
           nrow = C, byrow = TRUE),
    matrix(c(-.02, .12, 0,
             -.60, .66, 0,
             -1.26, -1.23,0),
           nrow = C, byrow = TRUE)
)

### List of betas for ILR model
beta_ILR_list <- list(
    matrix(c(.11, -.09,
             .99, .12,
             -.04, .60),
           nrow = C, byrow = TRUE),
    matrix(c(.16,-.04,
             .96, -0.03,
             0.03, 1.01),
           nrow = C, byrow = TRUE)
)


### Functions for different data generating mechanisms
dir_data_gen <- function(phi, mu) {
    ### mu is a N x C matrix with ith row specifying mean of observation i
    ### row sums must be 1
    N <- nrow(mu)
    C <- ncol(mu)
    yout <- matrix(NA, nrow = N, ncol = C)
    for(i in 1:N) {
        yout[i,] <- gtools::rdirichlet(1, alpha = phi * mu[i,])
    }
    return(yout)
}

mult_data_gen <- function(nmax, mu) {
    ### mu is a N x C matrix with ith row specifying mean of observation i
    ### row sums must be 1
    N <- nrow(mu)
    C <- ncol(mu)
    yout <- matrix(NA, nrow = N, ncol = C)
    n <- sample(1:nmax, N, replace = TRUE)
    for(i in 1:N) {
        mult_samples <- rmultinom(1, n[i], mu[i,])
        yout[i,] <- mult_samples[,1] / sum(mult_samples[,1])
    }
    return(yout)
}

dir_mult_data_gen <- function(nmax, phi, mu) {
    N <- nrow(mu)
    C <- ncol(mu)
    yout <- matrix(NA, nrow = N, ncol = C)
    n <- sample(1:nmax, N, replace = TRUE)
    for(i in 1:N) {
        pi <- gtools::rdirichlet(1, phi * mu[i,])
        mult_samples <- rmultinom(1, n[i], pi)
        yout[i,] <- mult_samples[,1] / sum(mult_samples[,1])
    }
    return(yout)
}

ilr_data_gen <- function(muIlr, sigma2) {
    N <- nrow(muIlr)
    K <- ncol(muIlr)
    yilr <- matrix(NA, nrow = N, ncol = K)
    for(k in 1:K) {
        yilr[,k] <- rnorm(N, muIlr[,k], sd = sqrt(sigma2)) 
    }
    y <- ilrInv(yilr)
    return(y)
}

### Function to remove zeros
removeZero <- function(compmat, eps = 1e-6) {
    compmat[compmat == 0] <- eps
    compmat <- compmat/rowSums(compmat)
    return(compmat)
}

all_settings <- expand.grid(N = c(100, 250, 500, 1000),
                            param_index = c(1,2),
                            true_model = c('compreg', 'ILR', 'mflr'),
                            datagen = c('dirichlet', 'multinomial', 'dirichlet-multinomial', 'ilr-normal'))
setting <- all_settings[index,]

### Generate covariates for evaluating predictions
set.seed(123)
Xtest <- gtools::rdirichlet(10000, rep(1, C))
XtestIlr <- ilr(Xtest)
XtestIlrIntercept <- cbind(1, XtestIlr)

### Generate true mu for evaluating predictions
if(setting$true_model == "ILR") {
    beta <- beta_ILR_list[[setting$param_index]]
    mu_ilr_true <- XtestIlrIntercept %*% beta
    mu_true <- ilrInv(mu_ilr_true)
} else if(setting$true_model  == "compreg"){
    mu_true <- Xtest %*% M_list[[setting$param_index]]
} else {
    beta <- beta_logit_list[[setting$param_index]]
    mu_true <- t(apply(XtestIlrIntercept %*% beta,1,softmax))
}

N <- setting$N
Nsims <- 10000
kld_df <- do.call(rbind, lapply(1:Nsims, function(sim) {
    ### Generate covariates
    X <- gtools::rdirichlet(N, rep(1, C))
    XIlr <- ilr(X)
    XIlrIntercept <- cbind(1, XIlr)
    
    ### Get mu for observed X
    if(setting$true_model == "ILR") {
        beta <- beta_ILR_list[[setting$param_index]]
        mu_ilr <- XIlrIntercept %*% beta
        mu <- ilrInv(mu_ilr)
    } else if(setting$true_model  == "compreg"){
        mu <- X %*% M_list[[setting$param_index]]
    } else {
        beta <- beta_logit_list[[setting$param_index]]
        mu <- t(apply(XIlrIntercept %*% beta,1,softmax))
    }
    
    ### Generate outcome
    if(setting$datagen == 'dirichlet') {
        yout <- dir_data_gen(10, mu)
    } else if(setting$datagen == 'multinomial') {
        yout <- mult_data_gen(30, mu)
    } else if(setting$datagen == 'dirichlet-multinomial') {
        yout <- dir_mult_data_gen(30, 10, mu)
    } else {
        yout <- ilr_data_gen(ilr(mu), 1)
    }
    ### CompReg
    compreg_est <- pairedCompReg(yout, X)
    compreg_pred <- Xtest %*% compreg_est
    compreg_kld <- kldMat(mu_true, compreg_pred)
    ### Ilr
    ilr_est <- lm(ilr(removeZero(yout)) ~ ilr(X))
    ilr_pred <- ilrInv(XtestIlrIntercept %*% ilr_est$coefficients)
    ilr_kld <- kldMat(mu_true, ilr_pred)
    ### MFLR
    mflr_est <- fixedUpdateMFLR(yout, XIlrIntercept)
    mflr_pred <- t(apply(XtestIlrIntercept %*% mflr_est,1,softmax))
    mflr_kld <- kldMat(mu_true, mflr_pred)
    out_df <- data.frame(kld = c(compreg_kld, ilr_kld, mflr_kld),
                         model = c("compreg", "ilr", "mflr"),
                         simnum = sim)
    return(out_df)
}))
results_df <- data.frame(kld_df, setting)

saveRDS(results_df, output_file)
quit('no')

    
    

