library(gtools)

softmax <- function(x){exp(x)/sum(exp(x))}

negll <- function(THETA, y, X)
{
    K <- ncol(y)
    negll <- 0
    PROBS <- t(apply(X %*% THETA,1,softmax))
    for(k in 1:K){
        negll <- negll - sum(y[,k]*log(PROBS[,k]) )
    }
    return(negll)
}

fixedUpdateMFLR <- function(y, X, delta = .0001) {
    d <- ncol(X) - 1
    K <- ncol(y)
    THETA.1 <- matrix(0,d+1,K)
    ETA.1 <- t(apply(X %*% THETA.1,1,softmax))
    ETA.0 <- 1 # let the loop begin
    iter <- 1
    negll.vals <- numeric()
    negll.vals[1] <- negll(THETA.1,y,X)
    # Build B inverse:
    B.inv <- - 2 * (diag(K-1) + matrix(1,K-1,K-1)) %x% solve( t(X) %*% X )
    # allow algorithm to begin
    SS <- 1
    while(max(abs(SS)) > delta)
    {
        ETA.0 <- ETA.1
        THETA.0 <- THETA.1
        ETA.1 <- t(apply(X %*% THETA.0,1,softmax))
        # Build score function SS
        SS <- numeric((K-1)*(d+1))
        for(k1 in 1:(K-1)) {
            ind1 <- ((k1-1)*(d+1) + 1):(k1*(d+1))
            SS[ind1] <- t(X) %*% (y[,k1]*ETA.1[,K]  - y[,K]*ETA.1[,k1])
        }
        # Perform update
        Theta.0.long <- as.vector(THETA.0[,-K])
        Theta.1.long <- Theta.0.long - B.inv %*% SS
        THETA.1 <- cbind(matrix(Theta.1.long,ncol=K-1),rep(0,d+1))
        ETA.1 <- t(apply(X %*% THETA.1,1,softmax))
        negll.vals[iter+1] <- negll(THETA.1,y,X)
        iter <- iter + 1
    }
    return(THETA.1)
}


