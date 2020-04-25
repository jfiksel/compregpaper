index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_independence_test_loss_function"
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

logLikILR <- function(y, ilrFit) {
    mu_ilr <- ilrFit$fitted.values[1:nrow(y),]
    sigma_ilr <- sigma(ilrFit)
    yIlr <- ilr(y)[1:nrow(y),]
    ll <- 0
    for(j in 1:ncol(yIlr)) {
        ll <- ll + sum(dnorm(yIlr[,j], mu_ilr[,j], sigma_ilr[j], log = TRUE))
    }
    return(ll)
}

removeZero <- function(compmat, eps = 1e-6) {
    compmat[compmat == 0] <- eps
    compmat <- compmat/rowSums(compmat)
    return(compmat)
}


permutationCompReg <- function(yout, ypred, nperms = 1000, init.seed = NULL) {
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
    y_avg <- colMeans(yout)
    ll_null <- sum(sweep(yout, MARGIN = 2, log(y_avg), `*`))
    
    llr_obs <- ll_full - ll_null
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- sample(1:nrow(ypred))
        ypred_perm <- ypred[perm_index,]
        m_perm <- pairedCompReg(yout, ypred_perm)
        ll_perm <- logLikComp(yout, ypred_perm, m_perm)
        permut_stats[i] <- ll_perm - ll_null
    }
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}

permutationILR <- function(yout, ypred, nperms = 1000, init.seed = NULL) {
    if(is.null(init.seed)){
        set.seed(123)
    } else {
        set.seed(init.seed)
    }
    ### Get observed beta
    yout <- removeZero(yout)
   # XIlr <- ilr(ypred)
   # XIlrIntercept <- cbind(1, XIlr)
    ilr_obs <- lm(ilr(yout) ~ ilr(ypred))
    llr_obs <- logLikILR(yout, ilr_obs)
    ### Get do permutation test
    permut_stats <- rep(NA, nperms)
    for(i in 1:nperms) {
        perm_index <- sample(1:nrow(ypred))
        ypred_perm <- ypred[perm_index,]
        ilr_perm <- lm(ilr(yout) ~ ilr(ypred_perm))
        ll_perm <- logLikILR(yout, ilr_perm)
        permut_stats[i] <- ll_perm 
    }
    pval <- mean(permut_stats >= llr_obs)
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

### Beta list for simulations
beta_list <- list(
    matrix(c(1, -2,
             2, -1,
             -1, 2),
           nrow = C, byrow = TRUE),
    matrix(c(1, -2,
             .333, -.333,
             -.333, .333),
           nrow = C, byrow = TRUE),
    matrix(c(1, -2,
             2, -1,
             0, 0),
           nrow = C, byrow = TRUE),
    matrix(c(1, -2,
             0, 0,
             0, 0),
           nrow = C, byrow = TRUE)
)

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
nperms <- 1000
#nperms <- 50
### First do using compositional regression model
pval_list_compreg <- lapply(c(100,250, 500, 1000), function(n) {
    pvals <- lapply(seq_along(M_list), function(M_index) {
        set.seed(seeds[index])
        ### Sample x from uniform dirichlet
        x <- gtools::rdirichlet(n, rep(1, C)) 
        M <- M_list[[M_index]]
        mu <- x %*% M
        ### Simulate y from dirichlet
        ydir <- dir_data_gen(10, mu)
        ### Simulate y from multinomial
        ymult <- mult_data_gen(30, mu)
        ### Simulate y from dirichlet-multinomial
        ydirmult <- dir_mult_data_gen(30, 10, mu)
        
        ### Get associated pvals for compreg
        pdircompreg <- permutationCompReg(ydir, x, nperms = nperms)
        pmultcompreg  <- permutationCompReg(ymult, x, nperms = nperms)
        pdirmultcompreg  <- permutationCompReg(ydirmult, x, nperms = nperms)
        p_compreg <- data.frame(p = c(pdircompreg, pmultcompreg, pdirmultcompreg),
                                datagen = c("dirichlet", "multinomial", "dirichlet-multinomial"),
                                Mindex = M_index,
                                n = n,
                                seed.index = index,
                                method = 'compreg',
                                true_model = 'compreg')
        ### Get associated pvals for MFLR
        pdirILR <- permutationILR(ydir, x, nperms = nperms)
        pmultILR  <- permutationILR(ymult, x, nperms = nperms)
        pdirmultILR  <- permutationILR(ydirmult, x, nperms = nperms)
        p_ILR <- data.frame(p = c(pdirILR, pmultILR, pdirmultILR),
                                datagen = c("dirichlet", "multinomial", "dirichlet-multinomial"),
                                Mindex = M_index,
                                n = n,
                                seed.index = index,
                                method = 'ILR',
                                true_model = 'compreg')
        return(rbind(p_compreg, p_ILR))
    })
    return(do.call(rbind, pvals))
})
pval_df_compreg <- do.call(rbind, pval_list_compreg)

##########################
### Now using logit model model
pval_list_ILR <- lapply(c(100,250, 500, 1000), function(n) {
    pvals <- lapply(seq_along(beta_list), function(M_index) {
        #print(paste0("n:", n, ", ", "M_index:", M_index))
        set.seed(seeds[index])
        ### Sample x from uniform dirichlet
        x <- gtools::rdirichlet(n, rep(1, C)) 
        xIlr <- ilr(x)
        XIlrIntercept <- cbind(1, xIlr)
        beta <- beta_list[[M_index]]
        muIlr <- XIlrIntercept %*% beta
        mu <- ilrInv(muIlr)
        ### Simulate y from ILR 
        yilr <- ilr_data_gen(ilr(mu), 1)
        yilr <- yilr[1:nrow(yilr),]
        
        ### Get associated pvals for compreg
        pilrcompreg <- permutationCompReg(yilr, x, nperms = nperms)
        
        p_compreg <- data.frame(p = pilrcompreg,
                                datagen = "ilr-normal",
                                Mindex = M_index,
                                n = n,
                                seed.index = index,
                                method = 'compreg',
                                true_model = 'ILR')
        ### Get associated pvals for MFLR
        pilrILR <- permutationILR(yilr, x, nperms = nperms)
        p_ILR <- data.frame(p = pilrILR,
                             datagen = "ilr-normal",
                             Mindex = M_index,
                             n = n,
                             seed.index = index,
                             method = 'ILR',
                             true_model = 'ILR')
        return(rbind(p_compreg, p_ILR))
    })
    return(do.call(rbind, pvals))
})
pval_df_ILR <- do.call(rbind, pval_list_ILR)
pval_df <- rbind(pval_df_compreg, pval_df_ILR)
saveRDS(pval_df, output_file)
quit('no')
