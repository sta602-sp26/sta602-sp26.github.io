# MVN with missing data
## From Hoff Ch. 7

library(tidyverse)
library(mvtnorm)
library(monomvn)
library(coda)
Y = read_csv("https://sta602-sp26.github.io/data/Pima.csv") %>%
  as.matrix() 
colnames(Y) = NULL

## prior parameters
n = nrow(Y); p = ncol(Y)
# prior on theta
mu0 = c(120, 64, 26, 26); sd0 = (mu0 / 2)
L0 = matrix(.1, p, p)
diag(L0) = 1 
L0 = L0 * outer(sd0, sd0) # \Lambda_0
# prior on Sigma
nu0 = p + 2; 
S0 = L0
###

### starting values
Sigma = S0
Y.full = Y
O = 1 * (!is.na(Y)) # indices for observe values of Y

for(j in 1:p) {
  Y.full[is.na(Y.full[,j]),j]<-mean(Y.full[,j],na.rm=TRUE)
} 

### Gibbs sampler
THETA <- SIGMA <- Y.MISS <- NULL
set.seed(360)

for(s in 1:1000) {
  
  ###update theta
  ybar <- apply(Y.full, 2 , mean)
  Ln <- solve(solve(L0) + n * solve(Sigma))
  mun <- Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% ybar)
  theta <- rmvnorm(1, mun, Ln)
  ###
  
  ###update Sigma
  Sn <- S0 + (t(Y.full) - c(theta)) %*% t(t(Y.full) - c(theta))
  Sigma <- rwish(nu0 + n, solve(Sn))
  ###
  
  ###update missing data
  for(i in 1:n) { 
    b <- (O[i, ] == 0)
    a <- (O[i, ] == 1)
    if( sum(b) != 0) {
      iSa <- solve(Sigma[a, a])
      beta.j <- Sigma[b, a] %*% iSa
      s2.j   <- Sigma[b, b] - Sigma[b, a] %*% iSa %*% Sigma[a, b]
      theta.j <- theta[b] + beta.j %*% (as.matrix(Y.full[i, a]) - theta[a])
      Y.full[i, b] <- rmvnorm(1, theta.j, s2.j)
    }
  }
  
  ### save results
  THETA<-rbind(THETA,theta)
  SIGMA<-rbind(SIGMA,c(Sigma))
  Y.MISS<-rbind(Y.MISS, Y.full[O==0] )
  ###
  
  if(s %% 250 == 0 | s == 1) {
    cat(s,theta,"\n")
  }
}

## Summaries

#### Posterior mean
apply(THETA,2,mean)

# posterior mean correlation matrix
COR <- array( dim=c(p,p,1000) )
for(s in 1:1000)
{
  Sig<-matrix( SIGMA[s,] ,nrow=p,ncol=p)
  COR[,,s] <- Sig/sqrt( outer( diag(Sig),diag(Sig) ) )
}

apply(COR,c(1,2),mean)

# effective sample size of THETA 
apply(THETA, 2, effectiveSize)

# other summaries
dim(Y.MISS)
Y.full[O==0]
colMeans(Y.MISS)
