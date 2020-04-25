index <- as.numeric(commandArgs(trailingOnly = TRUE))
output_dir <- "cluster_output_independence_test_multinomial"
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
output_file <- file.path(output_dir, paste0("run-", index, ".rds"))
if(file.exists(output_file)){
    quit('no')
}

source(file.path("..", "scripts", "compRegression.R"))

#### Function to do permutation test

logLikComp <- function(y, x, M) {
    mu <- x %*% M
    return(sum(y*log(mu)))
}


permutationCompReg <- function(yout, ypred, nperms = 1000, init.seed = NULL) {
    if(is.null(init.seed)){
        set.seed(123)
    } else {
        set.seed(init.seed)
    }
    ### Get observed M
    m_obs <- pairedCompReg(yout, ypred, accelerate = FALSE)
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
        m_perm <- pairedCompReg(yout, ypred_perm, accelerate = FALSE)
        ll_perm <- logLikComp(yout, ypred_perm, m_perm)
        permut_stats[i] <- ll_perm - ll_null
    }
    pval <- mean(permut_stats >= llr_obs)
    return(pval)
}

### set seed for simulations
set.seed(123)
seeds <- sample(-1e6:1e6, 10000, replace = F)
C <- 3
#########################
pvals <- lapply(c(100, 250, 500, 1000), function(n) {
    set.seed(seeds[index])
    ### Sample x from uniform dirichlet
    x <- t(rmultinom(n, 1, rep(1, C)))
    ### Simulate y from multinomial
    y <- t(rmultinom(n, 1, rep(1, C)))
    
    ### Get associated pvals
    p <- permutationCompReg(y, x)
    return(data.frame(p = p,
                      n = n,
                      seed.index = index))
})
pval_df <- do.call(rbind, pvals)
saveRDS(pval_df, output_file)
quit('no')
