### set seed for simulations
set.seed(123)
seeds <- sample(-1e6:1e6, 10000, replace = F)
C <- 3
#########################
p_list <- lapply(1:10000, function(index) {
    pvals <- lapply(c(100, 250, 500, 1000), function(n) {
        set.seed(seeds[index])
        ### Sample x from uniform dirichlet
        x <- t(rmultinom(n, 1, rep(1, C)))
        ### Simulate y from multinomial
        y <- t(rmultinom(n, 1, rep(1, C)))
        xcat <- max.col(x)
        ycat <- max.col(y)
        
        ### Get associated pvals
        p <- chisq.test(table(xcat, ycat))$p.value
        return(data.frame(p = p,
                          n = n,
                          seed.index = index))
    })
    return(do.call(rbind, pvals))
})
p_df <- do.call(rbind, p_list)

for(i in c(1000, 2000, 5000, 10000)) {
    print(p_df %>% filter(seed.index <= i) %>% group_by(n) %>% summarize(reject = mean(p <= .05)) )
}


p_list <- lapply(1:10000, function(index) {
    pvals <- lapply(c(100, 250, 500, 1000), function(n) {
        set.seed(seeds[index])
        ### Sample x from uniform dirichlet
        x <- t(rmultinom(n, 1, c(.45, .2, .35)))
        ### Simulate y from multinomial
        y <- t(rmultinom(n, 1, c(.25, .6, .15)))
        xcat <- max.col(x)
        ycat <- max.col(y)
        
        ### Get associated pvals
        p <- chisq.test(table(xcat, ycat))$p.value
        return(data.frame(p = p,
                          n = n,
                          seed.index = index))
    })
    return(do.call(rbind, pvals))
})
p_df <- do.call(rbind, p_list)

for(i in c(1000, 2000, 5000, 10000)) {
    print(p_df %>% filter(seed.index <= i) %>% group_by(n) %>% summarize(reject = mean(p <= .05)) )
}


