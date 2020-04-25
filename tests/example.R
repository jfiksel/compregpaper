library(gtools)
library(here)
source("compRegression.R")
#source(here("scripts", "compRegression.R"))

set.seed(123)
Dpred <- 3
Dout <- 2
N <- 300
### random ypred
ypred <- rdirichlet(N, rep(1, Dpred))
### Random coefficients 
M <- matrix(NA, nrow = Dpred, ncol = Dout)
for(i in 1:Dpred) {
    M[i,] <- rdirichlet(1, rep(1, Dout))
}

### Sample outcomes
yout <- matrix(NA, nrow = N, ncol = Dout)
for(i in 1:N) {
    mu <- 20*as.vector(t(M) %*% ypred[i,])
    yout[i,] <- rdirichlet(1, mu)
}

### run regression to get estimate of M
m_est <- pairedCompReg(yout, ypred)

### get predicted values

yhat <- matrix(NA, nrow = N, ncol = Dout)
for(i in 1:N) {
    yhat[i,] <- as.vector(t(M) %*% ypred[i,])
}


