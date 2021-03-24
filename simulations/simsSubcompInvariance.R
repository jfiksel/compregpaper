index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_compare_preds_revision"
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
M <- matrix(c(.65, .35, 0, 0, 0,
              0, .35, .65, 0, 0,
              .1, .1, .6, .1, .1,
              0, 0, 0, .8, .2,
              0, .4, 0, 0, .6),
            nrow = 5, byrow = TRUE) 
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
set.seed(123)
seeds <- sample(-1e6:1e6, size = 10000, replace = F)
est_list <- lapply(1:10000, function(i) {
    set.seed(seeds[i])
    X <- gtools::rdirichlet(250, rep(1, 5))
    mu <- X %*% M
    yout <- dir_data_gen(10, mu)
})
