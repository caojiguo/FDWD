
library(classiFunc)
data(DTI)
#citation, from classiFunc, acknowledgement

f <- function(x) any(is.na(x))
sum(apply(DTI$cca,1,f))
#================================
library(kerndwd)
library(fda)
library(fdapace)
source("FunctionalClassification.R")

n <- nrow(DTI$cca)
g <- ifelse(DTI$case==0,1,2)
table(g)

#windows()
pdf("DTIprofile.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
matplot(t(DTI$cca), type="l", col=g, lty=g, xlab="Distance Along Tract", ylab="Fractional Anisotropy")
dev.off()


#split 500 times, 10-fold CV 
B <- 500
classf.err <- matrix(NA, nrow=B, ncol=5)
for(b in 1:B)
{ #b=1
  set.seed(b)
  print(b)
  ID <- sample(c(rep(1:10,each=38), sample(1:10,size=2, prob=rep(.1,10))), size=n)
  pred.class <- matrix(NA,nrow=n, ncol=5)
  for(i in 1:10)
  { #i=1
    classf.DTI <- try(functional.classification(curves=t(DTI$cca[ID!=i,1:93]),groups=g[ID!=i],n1=sum(g[ID!=i]==1),
      n2=sum(g[ID!=i]==2),new.curves=t(DTI$cca[ID==i,1:93]),knmax=15,nbasis=40,norder=5,nfinegrid=301,tini=1,tend=93), T)
    if(inherits(classf.DTI, "try-error")) {next} else
    { 
      pred.class[ID==i,1:5] <- t(classf.DTI$groups.class)
    }
  }
  classf.err[b,1] <- mean(pred.class[,1]!=g, na.rm=T) #centroid
  classf.err[b,2] <- mean(pred.class[,2]!=g, na.rm=T) #functional quadratic discriminant analysis
  classf.err[b,3] <- mean(pred.class[,3]!=g, na.rm=T) #functional generalized linear model
  classf.err[b,4] <- mean(pred.class[,4]!=g, na.rm=T) # functional distance weighted discrimination with Gaussian kernel
  classf.err[b,5] <- mean(pred.class[,5]!=g, na.rm=T)  # functional distance weighted discrimination with exponential kernel
}

# save(classf.err, file="DTI.RData")
# load("DTI.RData")
colnames(classf.err) <- c("Cent","Quad","Logist","DWD_g","DWD_e")
#windows()
#pdf("DTImis.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
boxplot(classf.err,use.cols=T)
dev.off()