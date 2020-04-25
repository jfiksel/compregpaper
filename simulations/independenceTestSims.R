index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_independence_test"
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, paste0("run-", index, ".rds"))
if(file.exists(output_file)){
    quit('no')
}

source(file.path("..", "scripts", "compRegression.R"))
library(gtools)

#### Function to do permutation test


permutationCompReg <- function(yout, ypred, nperms = 1000, init.seed = NULL) {
    if(is.null(init.seed)){
        set.seed(123)
    } else {
        set.seed(init.seed)
    }
    ### Get observed M
    m_obs <- pairedCompReg(yout, ypred)
    ### Get observed test statistics
    y_avg <- colMeans(yout)
    chisq_obs <- sum(((m_obs - y_avg)^2)/y_avg)
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- sample(1:nrow(ypred))
        ypred_perm <- ypred[perm_index,]
        m_perm <- pairedCompReg(yout, ypred_perm)
        permut_stats[i] <- sum(((m_perm - y_avg)^2)/y_avg)
    }
    pval <- mean(permut_stats >= chisq_obs)
    return(pval)
}

### M list for simulations
### First is strong association
C <- 3
M1 <- .85*diag(C) + .05
### Second choice is moderate association
M2 <- matrix(c(.4, .3, .3,
               .3, .4, .3,
               .3, .3, .4),
             nrow = C, byrow = TRUE)
### Third choice of M represents association with first part of predictor
M3 <- matrix(c(.90, .05, .05,
               1/3, 1/3, 1/3,
               1/3, 1/3, 1/3),
             nrow = C, byrow = TRUE)
### Last choice is independence
M4 <- matrix(rep(1/C, C^2), nrow = C)
M_list <- list(M1, M2, M3, M4)

### Function to simulate from dirichlet
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

### Function to simulate from multinomial
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


setting.df <- expand.grid(n=c(100, 250, 500, 1000),
                          seed.index = 1:1000)
setting <- setting.df[index,]
### set seed for simulations
set.seed(123)
seeds <- sample(-1e6:1e6, 1000, replace = F)

#########################
pvals <- lapply(seq_along(M_list), function(M_index) {
    set.seed(seeds[setting$seed.index])
    ### Sample x from uniform dirichlet
    x <- rdirichlet(setting$n, rep(1, C)) 
    M <- M_list[[M_index]]
    mu <- x %*% M
    ### Simulate y from dirichlet
    ydir <- dir_data_gen(10, mu)
    ### Simulate y from multinomial
    ymult <- mult_data_gen(30, mu)
    ### Simulate y from dirichlet-multinomial
    ydirmult <- dir_mult_data_gen(30, 10, mu)
    
    ### Get associated pvals
    pdir <- permutationCompReg(ydir, x)
    pmult <- permutationCompReg(ymult, x)
    pdirmult <- permutationCompReg(ydirmult, x)
    return(data.frame(p = c(pdir, pmult, pdirmult),
                      datagen = c("dirichlet", "multinomial", "dirichlet-multinomial"),
                      Mindex = M_index,
                      n = setting$n,
                      seed.index = setting$seed.index))
})
pval_df <- do.call(rbind, pvals)
saveRDS(pval_df, output_file)
quit('no')
