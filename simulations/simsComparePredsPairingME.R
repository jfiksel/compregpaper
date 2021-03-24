index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_compare_preds_pairing_me"
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
M1 <- matrix(c(.7, .2, .1,
               0, .8, .2,
               .3, 0, .7),
             nrow = C, byrow = TRUE)

M2 <- matrix(c(.7, .2, .1,
               .4, .4, .2,
               .2, .2, .6),
             nrow = C, byrow = TRUE)

M_list <- list(M1, M2)

paired_measurements_data_gen <- function(M, N, C) {
    X <- Y <- matrix(0, nrow = N, ncol = C)
    prob_x <- gtools::rdirichlet(N, rep(1, C))
    for(i in 1:N) {
        r <- sample(1:50, 1)
        X[i,] <- rmultinom(1, size = r, prob = prob_x[i,])
        for(j in 1:C) {
            Y[i,] <- Y[i,] + rmultinom(1, size = X[i,j], prob = M[j,])
        }
        X[i,] <- X[i,] / sum(X[i,])
        Y[i,] <- Y[i,] / sum(Y[i,])
    }
    return(list(x = X, y = Y))
}

measurement_error_data_gen <- function(M, N, C) {
    X <- gtools::rdirichlet(N, rep(1, C))
    Y <- matrix(0, nrow = N, ncol = C)
    R <- 4000
    for(i in 1:N) {
        Z <- sample(1:C, 1, prob = X[i,])
        Y[i,] <- rmultinom(1, size = R, prob = M[Z,])
        Y[i,] <- Y[i,] / sum(Y[i,])
    }
    return(list(x = X, y = Y))
}

### Function to remove zeros
removeZero <- function(compmat, eps = 1e-6) {
    compmat[compmat == 0] <- eps
    compmat <- compmat/rowSums(compmat)
    return(compmat)
}

all_settings <- expand.grid(N = c(100, 250, 500, 1000),
                            param_index = c(1,2),
                            datagen = c('paired_aggregation', 'measurement_error'))
start_indices <- seq(1, 10000, by = 100) 
end_indices <- seq(100, 10000, by = 100)
indices <- data.frame(start = start_indices, end = end_indices )
all_settings <- do.call(rbind, lapply(1:nrow(all_settings), function(i) {
    setting <- all_settings[i,]
    return(cbind(setting, indices))
}))
setting <- all_settings[index,]

### generate seeds
set.seed(123)
my_seeds <- sample(-1e6:1e6, size = 10000, replace = FALSE)
### Generate data for evaluating predictions
X <- Y <- matrix(NA, nrow = 10000, ncol = C)
M_true <- M_list[[setting$param_index]]

if(setting$datagen == "paired_aggregation") {
    set.seed(123)
    paired_data <- paired_measurements_data_gen(M_true, 10000, C)
    Xtest <- paired_data$x
} else {
    set.seed(123)
    paired_data <- measurement_error_data_gen(M_true, 10000, C)
    Xtest <- paired_data$x
}
mu_true <- Xtest %*% M_list[[setting$param_index]]
XtestIlr <- ilr(removeZero(Xtest))
XtestIlrIntercept <- cbind(1, XtestIlr)

XtestALR <- alr(removeZero(Xtest))
XtestALRIntercept <- cbind(1, XtestALR)

N <- setting$N
start <- setting$start
end <- setting$end
kld_df <- do.call(rbind, lapply(start:end, function(sim) {
    set.seed(my_seeds[sim])
    if(setting$datagen == "paired_aggregation") {
        paired_data <- paired_measurements_data_gen(M_true, N, C)
        X <- paired_data$x
        yout <- paired_data$y
    } else {
        paired_data <- measurement_error_data_gen(M_true, N, C)
        X <- paired_data$x
        yout <- paired_data$y
    }
    XIlr <- ilr(removeZero(X))
    XIlrIntercept <- cbind(1, XIlr)
    
    XAlr <- alr(removeZero(X))
    XAlrIntercept <- cbind(1, XAlr)
    ### CompReg
    compreg_est <- pairedCompReg(yout, X)
    compreg_pred <- Xtest %*% compreg_est
    compreg_kld <- kldMat(mu_true, compreg_pred)
    ### Ilr
    ilr_est <- lm(ilr(removeZero(yout)) ~ ilr(removeZero(X)))
    ilr_pred <- ilrInv(XtestIlrIntercept %*% ilr_est$coefficients)
    ilr_kld <- kldMat(mu_true, ilr_pred)
    ### ALR
    alr_est <- lm(alr(removeZero(yout)) ~ alr(removeZero(X)))
    alr_pred <- alrInv(XtestALRIntercept %*% alr_est$coefficients)
    alr_kld <- kldMat(mu_true, alr_pred)
    ### MFLR
    mflr_est <- fixedUpdateMFLR(yout, XIlrIntercept)
    mflr_pred <- t(apply(XtestIlrIntercept %*% mflr_est,1,softmax))
    mflr_kld <- kldMat(mu_true, mflr_pred)
    out_df <- data.frame(kld = c(compreg_kld, ilr_kld, alr_kld, mflr_kld),
                         model = c("compreg", "ilr", "alr","mflr"),
                         simnum = sim)
    return(out_df)
}))
results_df <- data.frame(kld_df, setting)
saveRDS(results_df, output_file)
quit('no')