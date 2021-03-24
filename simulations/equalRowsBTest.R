index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_equal_rows_B_test"
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, paste0("run-", index, ".rds"))
if(file.exists(output_file)){
    quit('no')
}
#source(file.path("..", "scripts", "multFracReg.R"))
source(file.path("..", "scripts", "compRegression.R"))
library(compositions)
library(gtools)


#### Function to do permutation test

logLikComp <- function(y, x, M) {
    mu <- x %*% M
    return(sum(y*log(mu)))
}



removeZero <- function(compmat, eps = 1e-6) {
    compmat[compmat == 0] <- eps
    compmat <- compmat/rowSums(compmat)
    return(compmat)
}


permutationCompReg <- function(yout, ypred, pair_indices, nperms = 1000, init.seed = NULL) {
    if(is.null(init.seed)){
        set.seed(123)
    } else {
        set.seed(init.seed)
    }
    ### Get observed M
    m_obs <- pairedCompReg(yout, ypred)
    ### Get observed log-likelihood
    ll_full <- logLikComp(yout, ypred, m_obs)
    ### Get observed log-likelihood under null
    ynullpred <- ypred
    ynullpred[,pair_indices[1]] <- rowSums(ynullpred[,pair_indices[1:2]])
    ynullpred <- ynullpred[,-pair_indices[2]]
    mnull <- pairedCompReg(yout, ynullpred)
    ll_null <- logLikComp(yout, ynullpred, mnull)
    
    llr_obs <- ll_full - ll_null
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- matrix(NA, nrow = nrow(ypred), ncol = 2)
        perm_index[,1] <- rbinom(nrow(ypred), size = 1, prob = .5)
        perm_index[,2] <- 1 - perm_index[,1]
        perm_index <- perm_index + 1
        perm_index[,1] <- pair_indices[perm_index[,1]]
        perm_index[,2] <- pair_indices[perm_index[,2]]
        ypred_perm <- ypred
        for(r in 1:nrow(perm_index)) {
            ypred_perm[r,pair_indices] <-  ypred_perm[r,perm_index[r,]]
        }
        m_perm <- pairedCompReg(yout, ypred_perm)
        ll_perm <- logLikComp(yout, ypred_perm, m_perm)
        permut_stats[i] <- ll_perm - ll_null
    }
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}
C <- 3

M <- matrix(c(.90, .05, .05,
               1/3, 1/3, 1/3,
               1/3, 1/3, 1/3),
             nrow = C, byrow = TRUE)


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
    #n <- rep(nmax, N)
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
    #n <- rep(nmax, N)
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


#setting.df <- expand.grid(n=c(100, 250, 500, 1000),
#                          seed.index = 1:2000)
#setting <- setting.df[index,]
### set seed for simulations
set.seed(123)
seeds <- sample(-1e6:1e6, 10000, replace = F)

#########################
nperms <- 500
#nperms <- 50
pair_index_matrix <- matrix(c(1,2,
                              1,3,
                              2,3),
                            nrow = 3, byrow = TRUE)
### First do using compositional regression model
pval_list_compreg <- lapply(c(250), function(n) {
    set.seed(seeds[index])
    ### Sample x from uniform dirichlet
    x <- gtools::rdirichlet(n, rep(1, 3)) 
    mu <- x %*% M
    ### Simulate y from dirichlet
    ydir <- dir_data_gen(10, mu)
    ### Simulate y from multinomial
    ymult <- mult_data_gen(30, mu)
    ### Simulate y from dirichlet-multinomial
    ydirmult <- dir_mult_data_gen(30, 10, mu)
    pvals <- lapply(1:nrow(pair_index_matrix), function(pair_index) {
        ### Get associated pvals for compreg
        pdircompreg <- permutationCompReg(ydir, x, pair_index_matrix[pair_index,], nperms = nperms)
        pmultcompreg  <- permutationCompReg(ymult, x, pair_index_matrix[pair_index,], nperms = nperms)
        pdirmultcompreg  <- permutationCompReg(ydirmult, x, pair_index_matrix[pair_index,], nperms = nperms)
        p_compreg <- data.frame(p = c(pdircompreg, pmultcompreg, pdirmultcompreg),
                                datagen = c("dirichlet", "multinomial", "dirichlet-multinomial"),
                                pairindex = pair_index,
                                n = n,
                                seed.index = index,
                                method = 'compreg',
                                true_model = 'compreg')
        return(p_compreg)
    })
    return(do.call(rbind, pvals))
})
pval_df <- do.call(rbind, pval_list_compreg)
saveRDS(pval_df, output_file)
quit('no')
