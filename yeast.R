yeast <- read.table("yeast.txt")
ind <- apply(yeast, 1, function(x) sum(!is.finite(x)))
yeast <- yeast[ind<2,]
#===============================================
library(kerndwd)
library(fda)
library(fdapace)
library(locpol)
source("FunctionalClassification.R")
n <- nrow(yeast)
g <- yeast[,1]+1

windows()
#pdf("Yeastprofile.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
matplot(x=seq(0,119,length.out=18),y=t(as.matrix(yeast[,2:19])), type="l", col=g, lty=g, 
                          xlab="Time (minutes)", ylab="Gene Expression Level")
dev.off()


#presmoothing using local linear method
deg <- 1
kernel <- EpaK
est <- function(bw, dat, x) return(locPolSmootherC(dat$x,dat$y, x, bw, deg, kernel)$beta0)
t <- seq(0, 119, length.out=18)
t.grid <- seq(0,119,length.out=40)
curve <- as.matrix(yeast[,2:19])
smooth.yeast <- matrix(NA, nrow=n, ncol=40)
for(i in 1:n)
{ #i=1
  index1 <- which(is.finite(curve[i,]))
  cvBwSel <- regCVBwSelC(t[index1], as.numeric(curve[i,index1]), deg, kernel, interval = c(10, 20))
  smooth.yeast[i,] <-  est(cvBwSel, data.frame(x=t[index1], y=as.numeric(curve[i,index1])), t.grid)
}

save(smooth.yeast, file="yeast.RData")
load("yeast.RData")

#specify type is sparse when using pace. 

B <- 500
classf.err <- matrix(NA, nrow=B, ncol=5)
for(b in 1:B)
{ #b=436, starts from 437.
  set.seed(b)
  print(b)
  ID <- sample(c(rep(1:10,each=9), sample(1:10,size=1, prob=rep(.1,10))), size=n)
  pred.class <- matrix(NA,nrow=n, ncol=5)
  for(i in 1:10)
  { #i=1
    classf.yeast <- try(functional.classification(curves=t(smooth.yeast[ID!=i,]),groups=g[ID!=i],n1=sum(g[ID!=i]==1),
    n2=sum(g[ID!=i]==2),new.curves=t(smooth.yeast[ID==i,]),knmax=10,nbasis=30,norder=4,nfinegrid=301,tini=0,tend=119), T)
    if(inherits(classf.yeast, "try-error")) {next} else
    { 
      pred.class[ID==i,1:5] <- t(classf.yeast$groups.class)
    }
  }
  classf.err[b,1] <- mean(pred.class[,1]!=g, na.rm=T) #centroid
  classf.err[b,2] <- mean(pred.class[,2]!=g, na.rm=T) #functional quadratic discriminant analysis
  classf.err[b,3] <- mean(pred.class[,3]!=g, na.rm=T) #functional generalized linear model
  classf.err[b,4] <- mean(pred.class[,4]!=g, na.rm=T) # functional distance weighted discrimination with Gaussian kernel
  classf.err[b,5] <- mean(pred.class[,5]!=g, na.rm=T)  # functional distance weighted discrimination with exponential kernel
}


#save(classf.err, file="yeast.RData")
#load("yeast.RData"); 
colnames(classf.err) <- c("Cent","Quad","Logist","DWD_g","DWD_e")
#windows()
pdf("Yeastmis.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
boxplot(classf.err,use.cols=T)
dev.off()


apply(classf.err, 2, summary)
