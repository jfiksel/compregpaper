index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_sims_no_covar"
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
#source(here("scripts", "compRegression.R"))
source(file.path("..", "scripts", "compRegression.R"))

### Different data generating mechanisms
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

### three different transition matrices
C <- 4
M1 <- diag(C)
M2 <- matrix(c(1, 0, 0, 0,
               .65, .35, 0, 0,
               0, 0, 0.5, 0.5,
               0, 0, 0, 1),
             nrow = C, byrow = TRUE)
M3 <- .6*diag(C) + .1
M_list <- list(M1, M2, M3)

### Different setting data frames
dirichlet_setting <- expand.grid(dist = "dirichlet",
                                 phi = seq(10, 50, by = 10),
                                 nmax = NA,
                                 N = c(100, 250, 500, 1000),
                                 Mind = c(1, 2, 3))
mult_setting <- expand.grid(dist = "multinomial",
                                 phi = NA,
                                 nmax = seq(11, 51, by = 10),
                                 N = c(100, 250, 500, 1000),
                                 Mind = c(1, 2, 3))
mult_dir_setting <- expand.grid(dist = "multinomial-dirichlet",
                                phi = seq(10, 50, by = 10),
                                nmax = seq(11, 51, by = 10),
                                N = c(100, 250, 500, 1000),
                                Mind = c(1, 2, 3))

all_settings <-
    dirichlet_setting %>% 
    bind_rows(mult_setting) %>%
    bind_rows(mult_dir_setting)
setting <- all_settings[index,]

#### Generate predictor compositional vector and mean
N <- setting$N
M <- M_list[[setting$Mind]]
nsims <- 1000
#nsims <- 100
set.seed(123)
m_est_list <- lapply(1:nsims, function(i) {
    ypred <- gtools::rdirichlet(N, rep(1, C))
    mu <- ypred %*% M
    if(setting$dist == "dirichlet") {
        yout <- dir_data_gen(setting$phi, mu)
    } else if(setting$dist == "multinomial") {
        yout <- mult_data_gen(setting$nmax, mu)
    } else {
        yout <- dir_mult_data_gen(setting$nmax, setting$phi, mu)
    }
    estM <- pairedCompReg(yout, ypred)
    return(estM)
})

results_list <- list(setting = setting, m_est_list = m_est_list)
saveRDS(results_list, output_file)
quit('no')

