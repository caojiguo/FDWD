
#extract two classes, most difficult to classify
#"aa" as the vowel in "dark", and "ao" as the first vowel in "water"
dat <- read.table("phoneme.txt",header=T,sep=",")
dat.in <- dat[(dat$g=="ao")|(dat$g=="aa"), 2:258]
attributes(dat.in)$names <- NULL
#head(dat.in)
#class(dat.in[,257])
dat.in[,257] <- ifelse(dat.in[,257]=="ao", 1, 2)
#table(dat.in[,257])
#ao 1022 and aa 695

#=============================
#functional classifiers
library(kerndwd)
library(fda)
library(fdapace)
source("FunctionalClassification.R")
n <- nrow(dat.in)
g <- as.vector(unlist(dat.in[,257]))

#randomly select 100 curves to plot
set.seed(123)
ID <- sample(1:n,50)
windows()
#pdf("Phonemeprofile.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
matplot(t(as.matrix(dat.in[ID,1:256])), type="l", col=g[ID], lty=g[ID], xlab="Frequency", ylab="Log(periodogram)")
dev.off()

B <- 500
classf.err <- matrix(0, nrow=B, ncol=5)
for(b in 1:B)
{ #b=1
  set.seed(b)
  ID <- sample(c(rep(1:10,each=171), sample(1:10,size=7, prob=rep(.1,10))), size=n)
  pred.class <- matrix(NA,nrow=n, ncol=5)
  for(i in 1:10)
  { #i=1
    classf.phoneme <- try(functional.classification(curves=t(as.matrix(dat.in[ID!=i,1:256])),groups=g[ID!=i],n1=sum(g[ID!=i]==1),
      n2=sum(g[ID!=i]==2),new.curves=t(as.matrix(dat.in[ID==i,1:256])),knmax=15,nbasis=40,norder=5,nfinegrid=501,tini=1,tend=256), T)
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

# save(classf.err, file="phoneme.RData")
# load("phoneme.RData")
# apply(classf.err, 2, summary)
colnames(classf.err) <- c("Cent","Quad","Logist","DWD_g","DWD_e")
#windows()
#pdf("Phonememis.pdf")
par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(3, 1, 0), las=0)
boxplot(classf.err,use.cols=T)
dev.off()

