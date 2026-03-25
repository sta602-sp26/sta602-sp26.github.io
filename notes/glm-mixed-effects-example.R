# From Hoff Ch 11.
load(url("https://sta602-sp26.github.io/data/tumorLocation.RData"))
Y<-tumorLocation
xs<-seq(.05,1,by=.05)
X<-cbind(rep(1,ncol(Y)),poly(xs,deg=4,raw=TRUE))

### log-density of the multivariate normal distribution
ldmvnorm<-function(X,mu,Sigma,iSigma=solve(Sigma),dSigma=det(Sigma)) 
{
  Y<-t( t(X)-mu)
  sum(diag(-.5*t(Y)%*%Y%*%iSigma))  -
    .5*(  prod(dim(X))*log(2*pi) +     dim(X)[1]*log(dSigma) )
  
}
###


### sample from the multivariate normal distribution
rmvnorm<-function(n,mu,Sigma)
{
  p<-length(mu)
  res<-matrix(0,nrow=n,ncol=p)
  if( n>0 & p>0 )
  {
    E<-matrix(rnorm(n*p),n,p)
    res<-t(  t(E%*%chol(Sigma)) +c(mu))
  }
  res
}
###


### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}
###


### tumor example
XY.tumor<-dget("https://www2.stat.duke.edu/~pdh10/FCBS/Inline/XY.tumor")
###



### MCMC  

## mvnorm log density
ldmvnorm<-function(X,mu,Sigma,iSigma=solve(Sigma),dSigma=det(Sigma))
{
  Y<-t( t(X)-mu)
  sum(diag(-.5*t(Y)%*%Y%*%iSigma))  -
    .5*(  prod(dim(X))*log(2*pi) +     dim(X)[1]*log(dSigma) )
}

m<-nrow(Y) ; p<-ncol(X)


## priors
BETA = NULL
for(j in 1:21) {
  m1 = glm(Y[j,] ~ 0+X, family="poisson")
  BETA = rbind(BETA, m1$coefficients)
}

mu0 = colMeans(BETA)
S0<- 21*cov(BETA) ; eta0<-p+2
iL0<-iSigma<-solve(S0)

## MCMC
THETA.post<<-SIGMA.post<-NULL ; set.seed(1)
for(s in 1:50000) 
{
  
  ##update theta
  Lm<-solve( iL0 +  m*iSigma )
  mum<-Lm%*%( iL0%*%mu0 + iSigma%*%apply(BETA,2,sum) )
  theta<-t(rmvnorm(1,mum,Lm))
  
  ##update Sigma
  mtheta<-matrix(theta,m,p,byrow=TRUE)
  iSigma<-rwish(1,eta0+m, 
                solve( S0+t(BETA-mtheta)%*%(BETA-mtheta)) )
  
  ##update beta
  Sigma<-solve(iSigma) ; dSigma<-det(Sigma)
  for(j in 1:m)
  {
    beta.p<-t(rmvnorm(1,BETA[j,], .5*Sigma))
    
    lr<-sum( dpois(Y[j,],exp(X%*%beta.p),log=TRUE ) -
               dpois(Y[j,],exp(X%*%BETA[j,]),log=TRUE ) ) +
      ldmvnorm( t(beta.p),theta,Sigma,
                iSigma=iSigma,dSigma=dSigma ) -
      ldmvnorm( t(BETA[j,]),theta,Sigma,
                iSigma=iSigma,dSigma=dSigma )
    
    if( log(runif(1))<lr ) { BETA[j,]<-beta.p }
  }
  
  ##store some output
  if(s%%10==0)
  {  
    cat(s,"\n") 
    THETA.post<-rbind(THETA.post,t(theta)) 
    SIGMA.post<-rbind(SIGMA.post,c(Sigma)) 
  }
  
}

colMeans(THETA.post)
matrix(colMeans(SIGMA.post), nrow = 5, ncol = 5)

apply(THETA.post, FUN = coda::effectiveSize, MARGIN = 2)
