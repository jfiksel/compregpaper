### Goal is to assess variance for the different data generating mechanisms
library(gtools)
library(matrixStats)
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


set.seed(123)
X <- gtools::rdirichlet(100000, rep(1, C))
var_list <- lapply(1:2, function(m_index)  {
    mu <- X %*% M_list[[m_index]]
    ydir <- dir_data_gen(10, mu)
    ymult <- mult_data_gen(30, mu)
    ydirmult <- dir_mult_data_gen(30, 10, mu)
    return(matrix(c(colVars(ydir), colVars(ymult), colVars(ydirmult)),
                  nrow = 3, byrow = TRUE))
})
