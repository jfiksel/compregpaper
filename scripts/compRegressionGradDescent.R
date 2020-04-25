library(optimx)
library(gtools)
library(SQUAREM)

softmax <- function(x){exp(x)/sum(exp(x))}

ml2 <- function(A, XX, W, P, C, delta = .0001) {
    THETA.1 <- matrix(0,P,C)
    ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
    ETA.0 <- 1 # let the loop begin
    iter <- 1
    # Build B inverse:
    B.inv <- - 2 * (diag(C-1) + matrix(1,C-1,C-1)) %x% solve( t(XX) %*% XX )
    # allow algorithm to begin
    SS <- 1
    while(max(abs(SS)) > delta)
    {
        ETA.0 <- ETA.1
        THETA.0 <- THETA.1
        ETA.1 <- t(apply(XX %*% THETA.0,1,softmax))
        # Build score function SS
        SS <- numeric((C-1)*(P))
        for(k1 in 1:(C-1))
        {
            ind1 <- ((k1-1)*(P) + 1):(k1*(P))
            #SS[ind1] <- t(XX) %*% (A[,k1] * W[,k1] - ETA.1[,k1])
            #SS[ind1] <- t(XX) %*% (W[,k1]*(A[,k1]  - ETA.1[,k1]))
            SS[ind1] <- t(XX) %*% (W[,k1] * A[,k1]*ETA.1[,C]  - W[,C] * A[,C]*ETA.1[,k1])
        }
        # Perform update
        Theta.0.long <- as.vector(THETA.0[,-C])
        Theta.1.long <- Theta.0.long - B.inv %*% SS
        THETA.1 <- cbind(matrix(Theta.1.long,ncol=C-1),rep(0,P))
        ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
        iter <- iter + 1
    }
    return(THETA.1)
}

fixedUpdate <- function(A, XX, THETA.0, W, P, C, B.inv) {

    # Build B inverse:
    #B.inv <- - 2 * (diag(C-1) + matrix(1,C-1,C-1)) %x% solve( t(XX) %*% XX )
    # allow algorithm to begin
    ETA.1 <- t(apply(XX %*% THETA.0,1,softmax))
    # Build score function SS
    SS <- numeric((C-1)*(P))
    for(k1 in 1:(C-1))
    {
        ind1 <- ((k1-1)*(P) + 1):(k1*(P))
        #SS[ind1] <- t(XX) %*% (A[,k1] * W[,k1] - ETA.1[,k1])
        #SS[ind1] <- t(XX) %*% (W[,k1]*(A[,k1]  - ETA.1[,k1]))
        SS[ind1] <- t(XX) %*% (W[,k1] * A[,k1]*ETA.1[,C]  - W[,C] * A[,C]*ETA.1[,k1])
    }
    # Perform update
    Theta.0.long <- as.vector(THETA.0[,-C])
    Theta.1.long <- Theta.0.long - B.inv %*% SS
    THETA.1 <- cbind(matrix(Theta.1.long,ncol=C-1),rep(0,P))
    return(THETA.1)
}



weightedMultiNomML <- function(A, XX, W, P, C, delta = .001) {
    # initialize parameter values
    THETA.1 <- matrix(0,P,C)
    ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
    ETA.0 <- 1 
    
    # set convergence criterion
    iter <- 0
    while( max(abs(ETA.0 - ETA.1 )) > delta) # stop when fitted probabilities cease to change.
    {
        ETA.0 <- ETA.1
        THETA.0 <- THETA.1
        for( k in 1:(C-1)) {
            THETA.0[,k] <- THETA.1[,k] + 1 # let the loop begin
            while( max(abs(THETA.1[,k] - THETA.0[,k])) > delta) {
                THETA.0[,k] <- THETA.1[,k]
                p.k <- t(apply(XX %*% THETA.0,1,softmax))[,k]
                w.k <- p.k*(1-p.k)
                # Z.k <- XX %*% THETA.0[,k] + ( (Y == k) - p.k) / w.k
                Z.k <- XX %*% THETA.0[,k] + ( A[,k] * W[,k] - p.k) / w.k
                W.k <- diag(w.k)
                XX.w.k <- sqrt(W.k) %*% XX
                Z.w.k <- sqrt(w.k) * Z.k
                theta.1.k <- solve(t(XX.w.k) %*% XX.w.k) %*% t(XX.w.k) %*% Z.w.k
                THETA.1[,k] <- theta.1.k
            }
        }
        ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
        iter <- iter + 1
    }
    return(THETA.1)
}

weightedMultiGradDesc <- function(A, XX, W, P, C,THETA.0) {
    THETA.1 <- THETA.0
    for( k in 1:(C-1)) {
        p.k <- t(apply(XX %*% THETA.0,1,softmax))[,k]
        w.k <- p.k*(1-p.k)
        Z.k <- XX %*% THETA.0[,k] + ( A[,k] * W[,k] - p.k) / w.k
        #Z.k <- XX %*% THETA.0[,k] + ( V[,k] * W[,k] - p.k) / w.k
        W.k <- diag(w.k)
        XX.w.k <- sqrt(W.k) %*% XX
        Z.w.k <- sqrt(w.k) * Z.k
        theta.1.k <- solve(t(XX.w.k) %*% XX.w.k) %*% t(XX.w.k) %*% Z.w.k
        THETA.1[,k] <- theta.1.k
    }
    return(THETA.1)
}

weightedMultiNomML2 <- function(A, XX, W, P, C, delta = .0001) {
    # initialize parameter values
    THETA.1 <- matrix(0,P,C)
    ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
    ETA.0 <- 1 
    
    # set convergence criterion
    iter <- 0
    while( max(abs(ETA.0 - ETA.1 )) > delta) # stop when fitted probabilities cease to change.
    {
       # ETA.0 <- ETA.1
        THETA.0 <- THETA.1
        for( k in 1:(C-1)) {
            THETA.0[,k] <- THETA.1[,k] + 1 # let the loop begin
            while( max(abs(THETA.1[,k] - THETA.0[,k])) > delta) {
                THETA.0[,k] <- THETA.1[,k]
                ETA.0 <- t(apply(XX %*% THETA.0,1,softmax))
                S.k <- t(XX) %*% (W[,k] * A[,k]*ETA.1[,C]  - W[,C] * A[,C]*ETA.1[,k])
                w.k <- -W[,k] * A[,k] * ETA.1[,k] + W[,C] * A[,C] * (ETA.1[,k]) * (1-ETA.1[,k]) 
                W.k <- diag(w.k)
                H <- t(XX) %*% W.k %*% XX
                theta.1.k <- solve(t(XX.w.k) %*% XX.w.k) %*% t(XX.w.k) %*% Z.w.k
                THETA.1[,k] <- theta.1.k
            }
        }
        ETA.1 <- t(apply(XX %*% THETA.1,1,softmax))
        iter <- iter + 1
    }
    return(THETA.1)
}


compreg.em <- function(pars, A, X, G, C, P,B.inv) {
    betas <- pars
    dim(betas) <- c(P,C,C)
    betanew <- array(NA, dim = c(P,C,C))
    ### ### Need to compute weights W_[r,i,j]
    weights_array <- array(NA, dim = c(nrow(A), C, C))
    ### W_{rij} = weights_array[r,i,j]
    for(r in 1:N) {
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            M[i,] <- t(apply(X[r,]%*% betas[,,i],1,softmax))
        }
        for(j in 1:C) {
            weights <- sweep(t(as.matrix(G[r,])), MARGIN=2, M[,j], `*`)
            weights <- weights/rowSums(weights)
            weights_array[r,,j] <- weights
        } 
    }
    
    for(i in 1:C) {
        #betanew[,,i] <- weightedMultiNomML(A, XX = X, W = weights_array[,i,], P, C)
       # betanew[,,i] <- ml2(A, XX = X, W = weights_array[,i,], P, C)
        betanew[,,i] <- fixedUpdate(A, XX = X, betas[,,i], W = weights_array[,i,], P, C, B.inv)
        #betanew[,,i] <- weightedMultiGradDesc(A, XX = X, W = weights_array[,i,], P, C, betas[,,i])
    }
    parsnew <- as.vector(betanew)
    return(parsnew)
}

loglik <- function(pars, A, G, X, C, P,B.inv) {
    ### First turn into array
    betas <- pars
    dim(betas) <- c(P,C, C)

    N <- nrow(A)
    ll <- 0
    for(r in 1:N) {
        ### Create individal M matrix
        M <- matrix(NA, nrow = C, ncol = C)
        for(i in 1:C) {
            M[i,] <- t(apply(X[r,]%*% betas[,,i],1,softmax))
        }
        mu <- t(M) %*% G[r,]
        ### calculate log likelihood
        for(j in 1:C) {
            ll <- ll + A[r,j] * log(mu[j])
        }
    }
    return(-ll)
}



set.seed(123)
C <- 2
N <- 1000
### Number of columns for X
P <- 2
#G <- rdirichlet(N, rep(1, C))
G <- t(rmultinom(N, 1, prob = rep(1/C, C)))

### for each i, beta_{i,j} has to be of length P

beta_array <- array(NA, dim = c(P, C, C))
for(i in 1:C) {
    beta_array[,,i] <- cbind(matrix(rnorm((P)*(C-1), sd = 2),P,C-1),rep(0,P))
}


X <- cbind(rep(1,N),matrix(rbinom(N, 1, .5)))

A <- matrix(NA, nrow = N, ncol = C)
for(r in 1:N) {
    indivM <- matrix(NA, nrow = C, ncol = C)
    for(i in 1:C) {
        indivM[i,] <- t(apply(X[r,]%*% beta_array[,,i],1,softmax))
    }
    A[r,] <- rdirichlet(1, 10 * (t(indivM) %*% G[r,])[,1])
}


par0 <- rep(0, prod(dim(beta_array)))
B.inv <- - 2 * (diag(C-1) + matrix(1,C-1,C-1)) %x% solve( t(XX) %*% XX )
test <- squarem(par = par0, fixptfn = compreg.em,
               objfn = loglik, control = list(tol = 1.e-8),
               A = A, X = X,G = G, C = C, P = P, B.inv = B.inv)


beta_est <- test$par
dim(beta_est) <- c(P, C, C)




set.seed(123)
N <-5000
P <- 2
C <- 3
XX <- cbind(rep(1,N),matrix(rnorm(n*(P-1)),N,P-1))
THETA <- cbind(matrix(rnorm((P)*(C-1)),P,C-1),rep(0,P))
ETA <- t(apply(XX %*% THETA,1,softmax))
Y <- apply(ETA,1,FUN=sample, size=1,x=1:C, replace=FALSE)
A <- matrix(NA, nrow = N, ncol = C)
for(r in 1:N) {
    for(i in 1:C) {
        if(Y[r]==i) {
            A[r,i] <-1
        } else {
            A[r,i]<-0
        }
    }
}

weightedMultiNomML(A,XX, W = matrix(1, N, C), P, C, delta = 1e-4)



set.seed(123)
N <-1000
P <- 2
C <- 4
XX <- cbind(rep(1,N),matrix(rnorm(N*(P-1)),N,P-1))
THETA <- cbind(matrix(rnorm((P)*(C-1)),P,C-1),rep(0,P))
ETA <- t(apply(XX %*% THETA,1,softmax))
Y <- apply(ETA,1,FUN=sample, size=1,x=1:C, replace=FALSE)
A <- matrix(0, nrow = N, ncol = C)
for(i in 1:N) {
    for(j in 1:C) {
        if(Y[i] == j) {
            A[i,j] <- 1
        }
    }
}

test1 <-weightedMultiNomML2(A,XX, W = matrix(1, N, C), P, C)
test2 <- ml2(A,XX, W = matrix(1, N, C), P, C)

