library(gtools)
library(SQUAREM)

compreg.em <- function(pars, V, G, C) {
    ### Note that pars should be a C^2 vector representing M
    M <- matrix(pars, ncol = C, nrow = C, byrow = FALSE)
    Mnew <- matrix(NA, nrow = C, ncol = C)
    ### First compute generic weights for E[Z_i | d = j]
    weights_array <- array(NA, dim = c(nrow(V), C, C))
    ### W_{rij} = weights_array[r,i,j]
    for(j in 1:C) {
        weights <- sweep(G, MARGIN=2, M[,j], `*`)
        weights <- weights/rowSums(weights)
        weights_array[,,j] <- weights
    }
    for(i in 1:C) {
        for(j in 1:C) {
            w_ij <- weights_array[,i,j]
            Mnew[i,j] <- sum(V[,j] * w_ij)
        }
    }
    Mnew[Mnew<=0] <- 1e-8
    Mnew <- Mnew/rowSums(Mnew)
    parsnew <- as.vector(Mnew)
    return(parsnew)
}

compreg.loglik <- function(pars, V, G, C) {
    N <- nrow(V)
    M <- matrix(pars, ncol = C, nrow = C, byrow = FALSE)
    #M[M==0] <- .0001
    #M <- M/rowSums(M)
    # if(sum(M<0) > 0) {
    #     print("cannot have negative entries")
    #     return(-Inf)
    # }
    mu <- G %*% M
    #mu[mu==0] <- .0001
    #mu <- mu/rowSums(mu)
    loglik <- 0
    for(r in 1:N) {
        for(j in 1:C) {
            loglik <- loglik + V[r,j] * log(mu[r,j] + 1)
        }
    }
    return(-loglik)
}


### Simulate data
set.seed(123)
C <- 5
N <- 1000
G <- rdirichlet(N, rep(1, C))
M_true <- matrix(c(.65, .35, 0, 0, 0,
                   0, .35, .65, 0, 0,
                   .1, .1, .6, .1, .1,
                   0, 0, 0, .8, .2,
                   0, .4, 0, 0, .6),
                 nrow = C, byrow = TRUE) 
M_true[M_true==0] <- .0001
M_true <- M_true/rowSums(M_true)
mu <- G %*% M_true
A <- matrix(NA, nrow = N, ncol = C)
for(i in 1:N) {
    A[i,] <- rdirichlet(1, 10 * mu[i,])
}

### Create V_matrix
round_preserve_sum <- function(x, digits = 2) {
    up <- 10 ^ digits
    x <- x * up
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    return(y / up)
}

A_round <- t(apply(A, 1, function(x) round_preserve_sum(x, 2)))
V <- A_round * 100


initM <- matrix(NA, nrow = C, ncol = C)
initM <- rdirichlet(C, rep(1, C))

compreg.maxlik <- squarem(par = as.vector(initM), fixptfn = compreg.em,
                           objfn = compreg.loglik, control = list(tol = 1.e-08),
                           V = V, G = G, C = C)
Mest <- matrix(compreg.maxlik$par, ncol = C, nrow = C, byrow = FALSE)
