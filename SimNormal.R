library(kerndwd)
library(fda)
library(fdapace)
library(mvtnorm)
library(locpol)
source("FunctionalClassification.R")

deg <- 1
kernel <- EpaK
est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg, kernel)$beta0)

#number of grid points (50, 200)
m <- 200
t.grid <- seq(0, 1, length.out=m)
t.finegrid <- seq(0,1,length.out=200)
mean.f <- function(x) x + sin(x)
mean.eval <- mean.f(t.grid)
eigenv <- exp(-(1:21)/2)
fourb <- create.fourier.basis(rangeval=range(t.grid), nbasis=21)
basismat <- eval.basis(t.grid, fourb)

#varanice of measurement error in X
sig.err <- 0.1  #variance 0.01
n.train <- 50  #training data with 50 or 100
n.test <- 500  #test data always 500
n <- n.train+n.test
B <- 500  #500 replicates
classf.err <- matrix(NA, nrow=B, ncol=5)

y <- numeric(n)
X <- matrix(0, nrow=n, ncol=200)
for(b in 1:B)
{ #b=2
  print(b)
  set.seed(b)
  xi <- rmvnorm(n, sigma=diag(eigenv))
  for(i in 1:n)
  { #i=1
    W <- mean.eval + basismat %*% xi[i,] + rnorm(m, sd=sig.err)
    cvBwSel <- regCVBwSelC(t.grid, W, deg, kernel, interval = c(0, 0.1))
    X[i,] <-  est(cvBwSel, data.frame(x=t.grid, y=W), t.finegrid)
    y[i] <-  sign(xi[i,1]^2+cos(xi[i,1]*xi[i,2])+exp(xi[i,2]^2)+log(abs(xi[i,3])+1)-3)
  }
  
  g <- 1/2*y+3/2
  classf.norm <- try(functional.classification(curves=t(as.matrix(X[1:n.train,])),groups=g[1:n.train],n1=sum(g[1:n.train]==1),
      n2=sum(g[1:n.train]==2),new.curves=t(as.matrix(X[(1+n.train):n,])),knmax=15,nbasis=40,norder=5,nfinegrid=401,tini=0,tend=1), T)
  if(inherits(classf.norm, "try-error")) {next} else
  { 
    pred.class <- classf.norm$groups.class
    classf.err[b,1] <- mean(pred.class[1,]!=g[(1+n.train):n]) #centroid
    classf.err[b,2] <- mean(pred.class[2,]!=g[(1+n.train):n]) #functional quadratic discriminant analysis
    classf.err[b,3] <- mean(pred.class[3,]!=g[(1+n.train):n]) #functional generalized linear model
    classf.err[b,4] <- mean(pred.class[4,]!=g[(1+n.train):n]) # functional distance weighted discrimination with Gaussian kernel
    classf.err[b,5] <- mean(pred.class[5,]!=g[(1+n.train):n]) # functional distance weighted discrimination with exponential kernel
    
  }
}
windows()
boxplot(classf.err, use.cols=T)
#save(err,file="NormP50N100.RData")
#load("NormP50N100.RData")

no.na <- function(x) {all(is.na(x))}
index <- apply(classf.err, 1, no.na)
err <- classf.err[!index,]
apply(err, 2, mean)
apply(err, 2, sd)

#repeat the above procedure for n.train=50 (p=50)
#save(classf.err,file="NormP50N50.RData")
#load("NormP50N50.RData")

#===========================================
#200 equally spaced time points
#repeat the above procedure for n.train=100 (p=200)
m <- 200
t.grid <- seq(0, 1, length.out=m)
mean.f <- function(x) x + sin(x)
mean.eval <- mean.f(t.grid)
eigenv <- exp(-(1:21)/2)
fourb <- create.fourier.basis(rangeval=range(t.grid), nbasis=21)
basismat <- eval.basis(t.grid, fourb)

n.train <- 50  #training data with 50 or 100
n.test <- 500  #test data always 500
n <- n.train+n.test
B <- 500  #500 replicates
classf.err <- matrix(NA, nrow=B, ncol=5)

y <- numeric(n)
X <- matrix(0, nrow=n, ncol=m)

for(b in 1:B)
{ #b=1
  print(b)
  set.seed(b)
  xi <- rmvnorm(n, sigma=diag(eigenv))
  for(i in 1:n)
  {
    X[i,] <- mean.eval + basismat %*% xi[i,] + rnorm(m, sd=sig.err)
    y[i] <-  sign(xi[i,1]^2+cos(xi[i,1]*xi[i,2])+exp(xi[i,2]^2)+log(abs(xi[i,3])+1)-3)
  }
  #table(y)
  g <- 1/2*y+3/2
  classf.norm <- try(functional.classification(curves=t(as.matrix(X[1:n.train,])),groups=g[1:n.train],n1=sum(g[1:n.train]==1),
    n2=sum(g[1:n.train]==2),new.curves=t(as.matrix(X[(1+n.train):n,])),knmax=15,nbasis=40,norder=5,nfinegrid=401,tini=0,tend=1), T)
  if(inherits(classf.norm, "try-error")) {next} else
  { 
    pred.class <- classf.norm$groups.class
    classf.err[b,1] <- mean(pred.class[1,]!=g[(1+n.train):n]) #centroid
    classf.err[b,2] <- mean(pred.class[2,]!=g[(1+n.train):n]) #functional quadratic discriminant analysis
    classf.err[b,3] <- mean(pred.class[3,]!=g[(1+n.train):n]) #functional generalized linear model
    classf.err[b,4] <- mean(pred.class[4,]!=g[(1+n.train):n]) # functional distance weighted discrimination with Gaussian kernel
    classf.err[b,5] <- mean(pred.class[5,]!=g[(1+n.train):n]) # functional distance weighted discrimination with exponential kernel
    
  }
}

windows()
boxplot(classf.err, use.cols=T)
#save(err,file="NormP200N100.RData")
#load("NormP200N100.RData")

#repeat the above for n.train=50
#save(err,file="NormP200N50.RData")
#load("NormP200N50.RData")