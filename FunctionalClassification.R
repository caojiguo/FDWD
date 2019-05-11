
###############################################################################################
###############################################################################################
# This script contains two main functions:
# 1. FunctionalMahalanobis is a function that computes the functional Mahalanobis semi-distance 
# between a set of functional observations and its sample mean for a given functional dataset. 
# This function is written to consider only the case of curves observed at the same time points, 
# although it is not difficult to modify this function to consider curves observed at different 
# time points. Also, the function is written to use b-splines basis to smooth the curves.
# 2. FunctionsClassification is a function that performs classification for functional samples 
# with all the methods described in the paper.
# Other auxiliarity functions that are called by FunctionsClassification are included in the 
# script.
###############################################################################################
###############################################################################################

###############################################################################################
# The function FunctionalMahalanobis computes the Mahalanobis distance between a set of functional 
# observations and its sample mean.
# The inputs are:
# 1. curves, an m x n matrix with n curves (columns) observed at m points (rows)
# 2. tini, initial point of the observation interval
# 3. tend, final point of the observation interval
# 4. k, is the number of principal components taken to compute the functional Mahalanobis semi-distance
# 5. nbasis, is the number of basis used to fit the discrete curves
# 6. norder, is the number of interior knots
# The outputs are:
# 1. funcMah, the functional Mahalanobis semi-distances for the set of curves.
# 2. sco, the set of functional principal components scores.
# 3. eigenval, the set of eigenvalues
###############################################################################################

FunctionalMahalanobis <- function(curves,k,tini=0,tend=1,nbasis=40,norder=6){
  
  # Initial values
  
  m <- dim(curves)[1] # Number of observed points
  tm <- seq(tini,tend,length.out=m) # Observation points in the interval [tini,tend]
  
  # Smoothing with B-splines
  
  basisbsplines <- create.bspline.basis(c(tini,tend),nbasis,norder) # Creation of the basis
  curves.smooth <- smooth.basis(tm,curves,basisbsplines) # Data smoothing
  curves.smooth.fd <- curves.smooth$fd # Object containing the relevant information of smoothing
  
  # Obtention of FPCs of the curves
  
  curves.smooth.fd.pcalist <- pca.fd(curves.smooth.fd,k) # Computation of the FPCs
  eigval <- curves.smooth.fd.pcalist$values # Eigenvalues of the covariance operator
  sco <- curves.smooth.fd.pcalist$scores # Scores
  
  # Compute the functional Mahalanobis semi-distance
  
  funcMah <- mahalanobis(sco,colMeans(sco),diag(eigval[1:k]))
  
  return(list("funcMah"=funcMah,"sco"=sco,"eigval"=eigval))
  
}

#############################################################################################
#############################################################################################
# The function functional.classification provides with the classification of a set of curves
# given another set of curves with group assigned.
# The inputs are:
# 1. curves, an m x n matrix with n curves (columns) observed at m points (rows)
# 2. groups, an 1 x n vector with the class of each function in curves (1 for group 1, 2 for group 2)
# 3. n1, is the number of curves in the first group
# 4. n2, is the number of curves in the second group
# 5. new.curves, an m x n.new vector with n.new curves (columns) observed at m points (rows) with no
# class assigned
# 6. kmax, is the maximum number of neighbours in the kNN algorithm
# 7. knmax, is the maximum number of pcs used with the PCA and MAH distances
# 8. nbasis, is the number of b-spline basis functions used to fit the discrete curves
# 9. norder, is the number of interior knots
# 10. nfinegrid, is the number of points used in the fine grid used to compute some of the inner points
# 11. tini, is the initial point of the interal
# 12. tend, is the final point of the interval

# The outputs are:
# 1. groups.class, an 19 x n.new matrix with the class assigned by the 19 for each curve in the set of new.curves
# The rows are sorted as follows:

#1.Centorid with Deglaile and Hall method.
#2.Functional quadratic Bayes classification rule.
#3.Functional glm.
#4.Functional DWD method, Gaussian kernel.   
#5.Functional DWD method, exponential kernel

# 2. k.FPCs, an 1 x 11 vector with the number of functional principal components used in the above methods
# The vector is sorted as follows:

# 1. Centorid with Deglaile and Hall method.
# 2. Functional quadratic Bayes classification rule.
# 3.Functional glm.
# 4.Functional DWD method, Gaussian kernel.   
# 5.Functional DWD method, exponential kernel
#############################################################################################
#############################################################################################

functional.classification <- function(curves,groups,n1,n2,new.curves,knmax,nbasis,norder,nfinegrid,tini,tend){
# curves=t(BeetleFly[ID!=i,1:512]);groups=g[ID!=i];n1=sum(g[ID!=i]==1)
# n2=sum(g[ID!=i]==2);new.curves=t(BeetleFly[ID==i,1:512]);knmax=10;nbasis=40;norder=5;nfinegrid=601;tini=1;tend=512
  
  m <- dim(new.curves)[1]
  n.new <- dim(new.curves)[2]
  
  # Technical objects  
  
  tm <- seq(tini,tend,length.out=m) # Observation points in [tini,tend]
  tmfinegrid <- seq(tini,tend,length.out=nfinegrid) # Observation points in the fine grid
  hgrid <- (tend - tini)/(nfinegrid - 1) # Length of the steps between observation points in the fine grid
 # ksNN <- seq(1,kmax,by=2) # Values of k for kNN
 # lksNN <- length(ksNN) # Number of values of k for kNN
  
  # B-spline smoothing of the two set of curves ignoring the group membership
  
  basisbsplines <- create.bspline.basis(c(tini,tend),nbasis,norder) # Creation of the basis
  
  curves.smooth <- smooth.basis(tm,curves,basisbsplines) # Data smoothing for curves
  curves.smooth.fd <- curves.smooth$fd # Object containing the relevant information for curves
  curves.smooth.grid <- eval.fd(tmfinegrid,curves.smooth.fd) # Smoothed curves evaluated at a fine grid
  
  new.curves.smooth <- smooth.basis(tm,new.curves,basisbsplines) # Data smoothing for new.curves
  new.curves.smooth.fd <- new.curves.smooth$fd # Object containing the relevant information for new.curves
  new.curves.smooth.grid <- eval.fd(tmfinegrid,new.curves.smooth.fd) # Smoothed new.curves evaluated at a fine grid
  
  # B-spline smoothing of the set of curves with group membership. Also, computation of the means of curves for each group
  
  curves.g1 <- curves[,groups==1] # Curves in group 1
  curves.g1.smooth <- smooth.basis(tm,curves.g1,basisbsplines) # Data smoothing for curves in group 1
  curves.g1.smooth.fd <- curves.g1.smooth$fd # Object containing the relevant information for curves in group 1
  curves.g1.smooth.grid <- eval.fd(tmfinegrid,curves.g1.smooth.fd) # Smoothed curves in group 1 evaluated at a fine grid
  mean.g1.fd <- mean(curves.g1.smooth.fd) # Mean of curves in group 1
  mean.g1.fd.grid <- eval.fd(tmfinegrid,mean.g1.fd) # Mean of curves in group 1 evaluated at the fine grid
  
  curves.g2 <- curves[,groups==2] # Curves in group 2
  curves.g2.smooth <- smooth.basis(tm,curves.g2,basisbsplines) # Data smoothing for curves in group 2
  curves.g2.smooth.fd <- curves.g2.smooth$fd # Object containing the relevant information for curves in group 2
  curves.g2.smooth.grid <- eval.fd(tmfinegrid,curves.g2.smooth.fd) # Smoothed curves in group 2 evaluated at a fine grid
  mean.g2.fd <- mean(curves.g2.smooth.fd) # Mean of curves in group 2
  mean.g2.fd.grid <- eval.fd(tmfinegrid,mean.g2.fd) # Mean of curves in group 2 evaluated at the fine grid
  
  # Obtention of FPCs of curves assuming the same covariance operator (The curves are sorted by groups). Also given
  # eigenvalues, scores of centered curves, eigenfunctions evaluated at the fine grid, scores of new.curves and
  # scores for the difference of the means in group 1 and 2
  
  curves.same <- cbind(curves.g1 - apply(curves.g1,1,mean),curves.g2 - apply(curves.g2,1,mean)) # Centered curves
  groups.same <- c(rep(1,n1),rep(2,n2)) # Groups sorted by group
  curves.same.smooth <- smooth.basis(tm,curves.same,basisbsplines) # Data smoothing for centered curves
  curves.same.smooth.fd <- curves.same.smooth$fd # Object containing the relevant information for centered curves
  curves.same.smooth.fd.pcalist <- pca.fd(curves.same.smooth.fd,knmax) # Computation of the FPCs
  
  eigval.curves.same <- curves.same.smooth.fd.pcalist$values # Eigenvalues of the common covariance operator
  eigfun.curves.same <- eval.fd(tmfinegrid,curves.same.smooth.fd.pcalist$harmonics) # Eigenfunctions at the fine grid
  scores.curves.same <- curves.same.smooth.fd.pcalist$scores # Scores
  
  scores.new.curves.same.g1 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g1.fd.grid)) %*% eigfun.curves.same) # Scores new.curves g1 assuming same covariance operator
  scores.new.curves.same.g2 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g2.fd.grid)) %*% eigfun.curves.same) # Scores new.curves g2 assuming same covariance operator
  
  dif.means.g1.g2 <- mean.g1.fd.grid - mean.g2.fd.grid # Difference of the means in group 1 and 2 evaluated at a fine grid
  scores.dif.means.same.g1.g2 <- as.vector(hgrid * (t(as.vector(dif.means.g1.g2)) %*% eigfun.curves.same)) # Scores of difference means
  
  # Obtention of FPCs of curves assuming a different covariance operator (The curves are sorted by groups). Also given
  # eigenvalues, scores of centered curves, eigenfunctions evaluated at the fine grid, scores of new.curves and
  # scores for the difference of the means in group 1 and 2
  
  groups.dif <- c(rep(1,n1),rep(2,n2)) # Groups sorted by group
  curves.g1.smooth.fd.pcalist <- pca.fd(curves.g1.smooth.fd,knmax) # Computation of the FPCs for curves in group 1
  eigval.curves.g1 <- curves.g1.smooth.fd.pcalist$values # Eigenvalues for curves in group 1
  scores.curves.g1 <- curves.g1.smooth.fd.pcalist$scores # Scores for curves in group 1
  eigfun.curves.g1 <- eval.fd(tmfinegrid,curves.g1.smooth.fd.pcalist$harmonics) # Eigenfunctions for curves in group 1 at the fine grid
  
  curves.g2.smooth.fd.pcalist <- pca.fd(curves.g2.smooth.fd,knmax) # Computation of the FPCs for curves in group 2
  eigval.curves.g2 <- curves.g2.smooth.fd.pcalist$values # Eigenvalues for curves in group 2
  scores.curves.g2 <- curves.g2.smooth.fd.pcalist$scores # Scores for curves in group 2
  eigfun.curves.g2 <- eval.fd(tmfinegrid,curves.g2.smooth.fd.pcalist$harmonics) # Eigenfunctions for curves in group 2 at the fine grid
  
  scores.new.curves.dif.g1 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g1.fd.grid)) %*% eigfun.curves.g1) # Scores new.curves g1 assuming different covariance operator
  scores.new.curves.dif.g2 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g2.fd.grid)) %*% eigfun.curves.g2) # Scores new.curves g2 assuming different covariance operator
  
  # Coefficients of basis expansions
  
  coef.1 <- as.matrix(coef(curves.g1.smooth))
  coef.2 <- as.matrix(coef(curves.g2.smooth))
  coef.new <- as.matrix(coef(new.curves.smooth))
  
  ##########################################################################################################################   
  # Computation of the classifications for new.curves
  ##########################################################################################################################   
  
  groups.class <- matrix(NA,nrow=5,ncol=n.new)
  k.FPCs <- matrix(NA,nrow=1,ncol=5)
  
  
  cenfpcmahhdsame <- cen.fpc.mah.hd.same(curves.smooth.grid,mean.g1.fd.grid,mean.g2.fd.grid,eigfun.curves.same,eigval.curves.same,
                                         scores.new.curves.same.g1,scores.new.curves.same.g2,scores.dif.means.same.g1.g2,
                                         groups.same,n1,n2,n.new,knmax,hgrid)
  groups.class[1,] <- cenfpcmahhdsame$groups.hd.same
  k.FPCs[1] <- cenfpcmahhdsame$k.hd.same
  
  
  
  Bayesdif <- Bayes.dif(curves.smooth.grid,mean.g1.fd.grid,mean.g2.fd.grid,eigfun.curves.g1,eigfun.curves.g2,
                        hgrid,scores.new.curves.same.g1,scores.new.curves.same.g2,groups,
                        eigval.curves.g1,eigval.curves.g2,n1,n2,n.new,knmax)
  
  groups.class[2,] <- Bayesdif$groups.Bayes.dif
  k.FPCs[2] <- Bayesdif$k.Bayes.dif
  
 # FPCA with PACE without considering grouping
  n <- n1+n2
  curves.train <- vector(mode="list", length=n); curves.test <- vector(mode="list", length=n.new)
  time.train <- vector(mode="list", length=n); time.test <- vector(mode="list", length=n.new)
  for(i in 1:n)
  {
    curves.train[[i]] <- curves[,i]
    time.train[[i]] <- tm
  }
  for(i in 1:n.new)
  {
    curves.test[[i]] <- new.curves[,i]
    time.test[[i]] <- tm
  }
  p <- list(dataType="Dense",plot=F,methodSelectK='FVE',FVEthreshold=0.999,numBins=NULL,verbose= T)
  res = FPCA(curves.train,time.train,p)
  K.max <- min(ncol(res$xiEst),20)
  xi_train <- res$xiEst[,1:K.max]
  xi_test= predict(res, newLy=curves.test, newLt=time.test, sigma2=res$sigma2, K=K.max)
  
  # curves.same <- cbind(curves.g1 - apply(curves,1,mean),curves.g2 - apply(curves,1,mean)) # Centered curves
  # #groups.same <- c(rep(1,n1),rep(2,n2)) # Groups sorted by group
  # curves.same.smooth <- smooth.basis(tm,curves.same,basisbsplines) # Data smoothing for centered curves
  # curves.same.smooth.fd <- curves.same.smooth$fd # Object containing the relevant information for centered curves
  # curves.same.smooth.fd.pcalist <- pca.fd(curves.same.smooth.fd,knmax) # Computation of the FPCs
  # 
  # eigval.curves.same <- curves.same.smooth.fd.pcalist$values # Eigenvalues of the common covariance operator
  # xi_train <- curves.same.smooth.fd.pcalist$scores # Scores
  # eigfun.curves.same <- eval.fd(tmfinegrid,curves.same.smooth.fd.pcalist$harmonics) # Eigenfunctions at the fine grid
  # 
  # # scores.new.curves.same.g1 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g1.fd.grid)) %*% eigfun.curves.same) # Scores new.curves g1 assuming same covariance operator
  # # scores.new.curves.same.g2 <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.g2.fd.grid)) %*% eigfun.curves.same) # Scores new.curves g2 assuming same covariance operator
  # # 
  # # scores.new.curves.same <- 
  # curves.smooth.fd <- curves.smooth$fd # Object containing the relevant information for curves
  # mean.curves.fd <- mean(curves.smooth.fd)
  # mean.curves.fd.grid <- eval.fd(tmfinegrid,mean.curves.fd) 
  # 
  # # mean.g1.fd <- mean(curves.g1.smooth.fd) # Mean of curves in group 1
  # # mean.g1.fd.grid <- eval.fd(tmfinegrid,mean.g1.fd) 
  # #curves.smooth.grid <- eval.fd(tmfinegrid,curves.smooth.fd)
  # xi_test <- hgrid * (t(new.curves.smooth.grid - as.vector(mean.curves.fd.grid)) %*% eigfun.curves.same)
 
  
  #functional glm by Hans Muller et al.
  fit.fglm <- tune.fglm(train.score=xi_train, train.group=groups-1, test.score=xi_test)
  
  groups.class[3,] <- fit.fglm$groups.fglm
  k.FPCs[3] <- fit.fglm$k.fglm
  
  
  #DWD with gaussian Kernel and exponential Kernel
  n <- n1+n2
  zeta.train <- matrix(0, nrow=n, ncol=ncol(xi_train))
  zeta.test <- matrix(0, nrow=n.new, ncol=ncol(xi_test))

  for(j in 1:K.max)
  { #j=1
    zeta.train[,j] <- pnorm(xi_train[,j], sd=sqrt(res$lambda[j]))
    zeta.test[,j] <- pnorm(xi_test[,j], sd=sqrt(res$lambda[j]))
  }
  #dim(zeta.train); dim(zeta.test)
  #q=c(0.01,0.1, 1, 2, 10)
  sigma <- 10^(seq(-3,1,length.out=5))
  #gaussin kernel
  #kern.gau <-  rbfdot
  fit.fdwd <- tune.fdwd(train.score=xi_train, train.group=2*groups-3, test.score=xi_test, kern=rbfdot, sig=sigma)
  groups.class[4,] <- fit.fdwd$groups.fdwd
  k.FPCs[4] <- fit.fdwd$k.fdwd

  #exponential kernel
  #kern.exp <-  laplacedot
  fit.fdwd <- tune.fdwd(train.score=xi_train, train.group=2*groups-3, test.score=xi_test, kern=laplacedot, sig=sigma)
  groups.class[5,] <- fit.fdwd$groups.fdwd
  k.FPCs[5] <- fit.fdwd$k.fdwd
  
  
  #head(groups.class)
  
  return(list("groups.class"=groups.class,"k.FPCs"=k.FPCs))
  
}
# err.centroid <- mean(groups.class[1,]!=g[ID==0])
# err.fqda <- mean(groups.class[2,]!=g[ID==0])
# err.fglm <- mean(groups.class[3,]!=g[ID==0])
# err.fdwd.gau <- mean(groups.class[4,]!=g[ID==0])
# err.fdwd.exp <- mean(groups.class[5,]!=g[ID==0])


#############################################################################################
# cen.fpc.mah.hd.same provides with the classification with centroid method with FPC and Mahalanobis
# with the same covariance operator, and with HD method
#############################################################################################

cen.fpc.mah.hd.same <- function(curves.smooth.grid,mean.g1.fd.grid,mean.g2.fd.grid,eigfun.curves.same,eigval.curves.same,
                                scores.new.curves.same.g1,scores.new.curves.same.g2,scores.dif.means.same.g1.g2,
                                groups.same,n1,n2,n.new,knmax,hgrid){
  
  n <- n1 + n2
  
  good.class.fpc <- matrix(NA,nrow=knmax,ncol=1)
  good.class.mah <- matrix(NA,nrow=knmax,ncol=1)
  good.class.hd <- matrix(NA,nrow=knmax,ncol=1)
  
  for (kn in 1:knmax){
    
    groups.fpc <- matrix(NA,nrow=n,ncol=1)
    groups.mah <- matrix(NA,nrow=n,ncol=1)
    groups.hd <- matrix(NA,nrow=n,ncol=1)
    
    for (i in 1 : n){
      
      scores.1 <- hgrid * (t(curves.smooth.grid[,i] - as.vector(mean.g1.fd.grid)) %*% eigfun.curves.same)
      scores.2 <- hgrid * (t(curves.smooth.grid[,i] - as.vector(mean.g2.fd.grid)) %*% eigfun.curves.same)
      
      dist.fpc.g1 <- sqrt(sum((scores.1[1:kn])^2))
      dist.fpc.g2 <- sqrt(sum((scores.2[1:kn])^2))
      if (dist.fpc.g2 > dist.fpc.g1){groups.fpc[i] <- 1}else{groups.fpc[i] <- 2}
      
      dist.mah.g1 <- sqrt(sum((scores.1[1:kn])^2/eigval.curves.same[1:kn]))
      dist.mah.g2 <- sqrt(sum((scores.2[1:kn])^2/eigval.curves.same[1:kn]))
      if (dist.mah.g2 > dist.mah.g1){groups.mah[i] <- 1}else{groups.mah[i] <- 2}
      
      dist.hd.g1 <- sqrt(sum(abs(scores.1[1:kn] * scores.dif.means.same.g1.g2[1:kn])/eigval.curves.same[1:kn]))
      dist.hd.g2 <- sqrt(sum(abs(scores.2[1:kn] * scores.dif.means.same.g1.g2[1:kn])/eigval.curves.same[1:kn]))
      if (dist.hd.g2 > dist.hd.g1){groups.hd[i] <- 1}else{groups.hd[i] <- 2}
      
    } 
    
    good.class.fpc[kn] <- sum(groups.fpc==groups.same)
    good.class.mah[kn] <- sum(groups.mah==groups.same)
    good.class.hd[kn] <- sum(groups.hd==groups.same)
    
  }
  
  k.fpc.same <- min(which(good.class.fpc==max(good.class.fpc)))
  k.mah.same <- min(which(good.class.mah==max(good.class.mah)))
  k.hd.same <- min(which(good.class.hd==max(good.class.hd)))
  
  groups.fpc.same <- matrix(NA,nrow=n.new,ncol=1)  
  groups.mah.same <- matrix(NA,nrow=n.new,ncol=1)
  groups.hd.same <- matrix(NA,nrow=n.new,ncol=1)
  
  for (i in 1 : n.new){
    
    dist.fpc.g1 <- sqrt(sum((scores.new.curves.same.g1[i,1:k.fpc.same])^2))
    dist.fpc.g2 <- sqrt(sum((scores.new.curves.same.g2[i,1:k.fpc.same])^2))
    if (dist.fpc.g2 > dist.fpc.g1){groups.fpc.same[i] <- 1}else{groups.fpc.same[i] <- 2}
    
    dist.mah.g1 <- sqrt(sum((scores.new.curves.same.g1[i,1:k.mah.same])^2/eigval.curves.same[1:k.mah.same]))
    dist.mah.g2 <- sqrt(sum((scores.new.curves.same.g2[i,1:k.mah.same])^2/eigval.curves.same[1:k.mah.same]))
    if (dist.mah.g2 > dist.mah.g1){groups.mah.same[i] <- 1}else{groups.mah.same[i] <- 2}
    
    dist.hd.g1 <- sqrt(sum(abs(scores.new.curves.same.g1[i,1:k.hd.same] * scores.dif.means.same.g1.g2[1:k.hd.same])/eigval.curves.same[1:k.hd.same]))
    dist.hd.g2 <- sqrt(sum(abs(scores.new.curves.same.g2[i,1:k.hd.same] * scores.dif.means.same.g1.g2[1:k.hd.same])/eigval.curves.same[1:k.hd.same]))
    if (dist.hd.g2 > dist.hd.g1){groups.hd.same[i] <- 1}else{groups.hd.same[i] <- 2}
    
  }
  
  return(list("groups.fpc.same"=groups.fpc.same,"groups.mah.same"=groups.mah.same,"groups.hd.same"=groups.hd.same,
              "k.fpc.same"=k.fpc.same,"k.mah.same"=k.mah.same,"k.hd.same"=k.hd.same))
  
}




#############################################################################################
# Bayes.dif provides with the classification with Bayes rule assuming a different covariance 
# operator
#############################################################################################

Bayes.dif <- function(curves.smooth.grid,mean.g1.fd.grid,mean.g2.fd.grid,eigfun.curves.g1,eigfun.curves.g2,
                      hgrid,scores.new.curves.same.g1,scores.new.curves.same.g2,groups,
                      eigval.curves.g1,eigval.curves.g2,n1,n2,n.new,knmax){
  
  n <- n1 + n2
  pi1 <- n1 / (n1 + n2)
  pi2 <- n2 / (n1 + n2)
  
  good.class <- matrix(NA,nrow=knmax,ncol=1)
  for (kn in 1:knmax){
    groups.l <- matrix(NA,nrow=n,ncol=1)
    for (i in 1 : n){
      scores.1 <- hgrid * (t(curves.smooth.grid[,i] - as.vector(mean.g1.fd.grid)) %*% eigfun.curves.g1)
      scores.2 <- hgrid * (t(curves.smooth.grid[,i] - as.vector(mean.g2.fd.grid)) %*% eigfun.curves.g2)
      distg1 <- sum((scores.1[1:kn])^2/eigval.curves.g1[1:kn]) + log(prod(eigval.curves.g1[1:kn])) - 2 * log(pi1)
      distg2 <- sum((scores.2[1:kn])^2/eigval.curves.g2[1:kn]) + log(prod(eigval.curves.g2[1:kn])) - 2 * log(pi2)
      if (distg2 > distg1){groups.l[i] <- 1}else{groups.l[i] <- 2}
    }
    good.class[kn] <- sum(groups.l==groups)
  }
  wkn <- min(which(good.class==max(good.class)))
  
  groups.Bayes.dif <- matrix(NA,nrow=n.new,ncol=1)
  for (i in 1 : n.new){
    distg1 <- sum((scores.new.curves.same.g1[i,1:wkn])^2/eigval.curves.g1[1:wkn]) + log(prod(eigval.curves.g1[1:wkn])) - 2 * log(pi1)
    distg2 <- sum((scores.new.curves.same.g2[i,1:wkn])^2/eigval.curves.g2[1:wkn]) + log(prod(eigval.curves.g2[1:wkn])) - 2 * log(pi2)
    if (distg2 > distg1){groups.Bayes.dif[i] <- 1}else{groups.Bayes.dif[i] <- 2}
  }
  
  return(list("groups.Bayes.dif"=groups.Bayes.dif,"k.Bayes.dif"=wkn))
  
}


#functional glm
tune.fglm <- function(train.score, train.group, test.score)
{ #train.score=xi_train; train.group=groups-1; test.score=xi_test
  d <- ncol(train.score)
  ntrain <- nrow(train.score)
  err <- matrix(NA, nrow=ntrain, ncol=d)
  ID <- sample(c(rep(1:5, each=floor(ntrain/5)), rep(1, ntrain-5*floor(ntrain/5))))
  for(j in 1:d)
  { 
    
    for(i in 1:5)
    { # j=2
      dat.glm <- data.frame(gr=train.group[ID!=i], x=train.score[ID!=i,1:j])
      fglm <- try(glm(gr~., data=dat.glm, family=binomial(link="logit")), T)
      if(inherits(fglm, "try-error")) {next} else
      {
        pred.1 <- ifelse(predict(fglm, newdata=data.frame(x=train.score[ID==i,1:j]), type="response") > 0.5, 1, 0)
        err[ID==i,j] <- mean(pred.1 != train.group[ID==i])
      }
      
    }  
  }
  d.opt <- which.min(apply(err, 2, mean, na.rm=T))
  dat.glm <- data.frame(gr=train.group, x=train.score[,1:d.opt])
  fglm <- glm(gr~., data=dat.glm,family=binomial(link="logit"))
  pred.group <- ifelse(predict(fglm, newdata=data.frame(x=test.score[,1:d.opt]), type="response") > 0.5, 2, 1)
  
  return(list("groups.fglm"=pred.group,"k.fglm"=d.opt))
}


#functional dwd
tune.fdwd <- function(train.score, train.group, test.score, kern, qvals, sig)
{
  #train.score=zeta.train; train.group=2*groups-3; test.score=zeta.test; kern=rbfdot; qvals=q; sig=sigma
  scale.train.score <- scale(train.score, center=T, scale=T)
  lambda = 10^(seq(3, -3, length.out=20))
  
  d <- ncol(train.score)
  L <- length(sig)
  err <- matrix(0, nrow=L, ncol=d)
  
 
 for(l in 1:L)
 { #l=1
   kfun=kern(sigma=sig[l])
    for(j in 2:d)
      { # j=2
       #print(c(l, j))
        m.cv <- cv.kerndwd(x=scale.train.score[,1:j], y=train.group, kern=kfun, qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
       # m.cv = tunedwd(x, y, kern, lambda=lambda, qvals=q, eps=1e-5, maxit=1e5)
       # err[l,j] <- m.cv$cvm.min
       err[l,j] <- m.cv$cvm.min
      }  
 }
 rowmin.index <- which.min(apply(err[,-1],1, min))
 sig.opt <- sig[rowmin.index]
 d.opt <- which.min(apply(err[,-1],2, min))
 
 kfun <- kern(sigma=sig.opt)
 m.cv <- cv.kerndwd(x=scale.train.score[,1:d.opt], y=train.group, kern=kfun, qval=1, lambda=lambda, eps=1e-5, maxit=1e5)
 
 foo <- kerndwd(x= scale.train.score[,1:d.opt],y=train.group, kern=kfun,qval=1,lambda=m.cv$lambda.min, eps=1e-5, maxit=1e5)
 scale.test.score <- scale(test.score, center=T, scale=T)
 pred.group <- predict(foo, kern=kfun, x=scale.train.score[,1:d.opt], newx=scale.test.score[,1:d.opt])
 pred.group <- (pred.group+3)/2
 return(list("groups.fdwd"=pred.group,"k.fdwd"=d.opt))
}