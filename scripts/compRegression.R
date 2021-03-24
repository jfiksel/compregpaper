library(optimx)
library(SQUAREM)
library(boot)

softmax <- function(x){exp(x)/sum(exp(x))}

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
        SS[ind1] <- t(XX) %*% (W[,k1] * A[,k1]*ETA.1[,C]  - W[,C] * A[,C]*ETA.1[,k1])
    }
    # Perform update
    if(P == 1) {
        Theta.0.long <- as.vector(THETA.0[-C]) 
    } else {
        Theta.0.long <- as.vector(THETA.0[,-C]) 
    }
    Theta.1.long <- Theta.0.long - B.inv %*% SS
    THETA.1 <- cbind(matrix(Theta.1.long,ncol=C-1),rep(0,P))
    return(THETA.1)
}

compreg.em <- function(pars, A, G, D1, D2) {
    ### Note that pars should be a C^2 vector representing M
    M <- matrix(pars, ncol = D1, nrow = D2, byrow = FALSE)
    Mnew <- matrix(NA, ncol = D1, nrow = D2)
    ### First compute generic weights for E[Z_i | d = j]
    weights_array <- array(NA, dim = c(nrow(A), D2, D1))
    ### W_{rij} = weights_array[r,i,j]
    for(j in 1:D1) {
        weights <- sweep(G, MARGIN=2, M[,j], `*`)
        weights <- weights/rowSums(weights)
        weights_array[,,j] <- weights
    }
    for(i in 1:D2) {
        for(j in 1:D1) {
            w_ij <- weights_array[,i,j]
            Mnew[i,j] <- sum(A[,j] * w_ij)
        }
    }
    Mnew[Mnew<=0] <- 1e-8
    Mnew <- Mnew/rowSums(Mnew)
    parsnew <- as.vector(Mnew)
    return(parsnew)
}


compreg.covar.em <- function(pars, A, X, G, C, P,B.inv) {
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
        betanew[,,i] <- fixedUpdate(A, XX = X, betas[,,i], W = weights_array[,i,], P, C, B.inv)
    }
    parsnew <- as.vector(betanew)
    return(parsnew)
}

compreg.loglik <- function(pars, A, G, D1, D2) {
    N <- nrow(A)
    M <- matrix(pars, ncol = D1, nrow = D2, byrow = FALSE)
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
        for(j in 1:D1) {
            loglik <- loglik + A[r,j] * log(mu[r,j] + .001)
        }
    }
    return(-loglik)
}

compreg.covar.loglik <- function(pars, A, G, X, C, P,B.inv) {
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

pairedCompReg <- function(yout, ypred, X = NULL, boot.ci = FALSE, R = 500, accelerate = FALSE) {
    #C <- ncol(yout)
    D1 <- ncol(yout)
    D2 <- ncol(ypred)
    if(is.null(X)) {
        par0 <- rep(1/D1, D1*D2)
        if(accelerate) {
            em_output <- squarem(par = par0, fixptfn = compreg.em,
                                 objfn = compreg.loglik, control = list(tol = 1.e-8),
                                 A = yout, G = ypred, D1 = D1, D2 = D2)  
        } else {
            em_output <- fpiter(par = par0, fixptfn = compreg.em,
                                objfn = compreg.loglik, control = list(tol = 1.e-8),
                                A = yout, G = ypred, D1 = D1, D2 = D2)
        }
        M_est <- em_output$par
        dim(M_est) <- c(D2, D1)
        return(M_est)
    } else{
        P <- ncol(X)
        # par0 <- rep(0, prod(dim(beta_array)))
        par0 <- rep(0, P*C^2)
        B.inv <- - 2 * (diag(C-1) + matrix(1,C-1,C-1)) %x% solve( t(X) %*% X )
        em_output <- squarem(par = par0, fixptfn = compreg.covar.em,
                             objfn = compreg.covar.loglik, control = list(tol = 1.e-8),
                             A = yout, X = X, G = ypred, C = C, P = P, B.inv = B.inv)
        beta_est <- em_output$par
        dim(beta_est) <- c(P, C, C)
        return(beta_est) 
    }
}


bootCompReg <- function(data, ypred, X = NULL, indices) {
    yout_boot <- data[indices,]
    ypred_boot <- ypred[indices,]
    D1 <- ncol(yout_boot)
    D2 <- ncol(ypred_boot)
    if(is.null(X)) {
        par0 <- rep(1/D1, D1*D2)
        em_output <- squarem(par = par0, fixptfn = compreg.em,
                             objfn = compreg.loglik, control = list(tol = 1.e-8),
                             A = yout_boot, G = ypred_boot, D1 = D1, D2 = D2)
        M_est <- em_output$par
        return(M_est)
    } else {
        P <- ncol(X)
        X_boot <- X[indices,]
        # par0 <- rep(0, prod(dim(beta_array)))
        par0 <- rep(0, P*C^2)
        B.inv <- - 2 * (diag(C-1) + matrix(1,C-1,C-1)) %x% solve( t(X) %*% X )
        em_output <- squarem(par = par0, fixptfn = compreg.covar.em,
                             objfn = compreg.covar.loglik, control = list(tol = 1.e-8),
                             A = yout_boot, X = X_boot, G = ypred_boot,
                             C = C, P = P, B.inv = B.inv)
        beta_est <- em_output$par
        return(beta_est) 
    }
}

