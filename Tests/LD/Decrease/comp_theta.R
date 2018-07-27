rm(list=ls())

# Fitness
setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo100/Fitness")
dataA1 <- read.table("tabA.csv",header=TRUE)
dataB1 <- read.table("tabB.csv",header=TRUE)
dataC1 <- read.table("tabC.csv",header=TRUE)
dataD1 <- read.table("tabD.csv",header=TRUE)

setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness")
dataA2 <- read.table("tabA.csv",header=TRUE)
dataB2 <- read.table("tabB.csv",header=TRUE)
dataC2 <- read.table("tabC.csv",header=TRUE)
dataD2 <- read.table("tabD.csv",header=TRUE)

par(mfrow=c(1,1))
plot(0:24,dataA1$mean_r,type='l',col="firebrick2",ylim=c(0,1.3),ylab="DL",xlab="t")
lines(0:24,dataA2$mean_r,type='l',col="firebrick2",lty=2)
lines(0:24,dataB1$mean_r,type='l',col="mediumorchid")
lines(0:24,dataB2$mean_r,type='l',col="mediumorchid",lty=2)
lines(0:24,dataC1$mean_r,type='l',col="dodgerblue4")
lines(0:24,dataC2$mean_r,type='l',col="dodgerblue4",lty=2)
lines(0:24,dataD1$mean_r,type='l',col="forestgreen")
lines(0:24,dataD2$mean_r,type='l',col="forestgreen",lty=2)
legend("topright",c("theta=0.5","theta=0.1","theta=0.01","theta=0","allof=100","allof=50"),
       lty=c(1,1,1,1,1,2),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1),cex=0.6)

rm(list=ls())

# Control
setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo100/Control")
dataA1 <- read.table("tabA.csv",header=TRUE)
dataB1 <- read.table("tabB.csv",header=TRUE)
dataC1 <- read.table("tabC.csv",header=TRUE)
dataD1 <- read.table("tabD.csv",header=TRUE)

setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Control")
dataA2 <- read.table("tabA.csv",header=TRUE)
dataB2 <- read.table("tabB.csv",header=TRUE)
dataC2 <- read.table("tabC.csv",header=TRUE)
dataD2 <- read.table("tabD.csv",header=TRUE)

plot(0:24,dataA1$mean_r,type='l',col="firebrick2",ylim=c(0,1.3),ylab="DL",xlab="t",
     main="Control, comparison allogamy rates")
lines(0:24,dataA2$mean_r,type='l',col="firebrick2",lty=2)
lines(0:24,dataB1$mean_r,type='l',col="mediumorchid")
lines(0:24,dataB2$mean_r,type='l',col="mediumorchid",lty=2)
lines(0:24,dataC1$mean_r,type='l',col="dodgerblue4")
lines(0:24,dataC2$mean_r,type='l',col="dodgerblue4",lty=2)
lines(0:24,dataD1$mean_r,type='l',col="forestgreen")
lines(0:24,dataD2$mean_r,type='l',col="forestgreen",lty=2)
legend("topright",c("theta=0.5","theta=0.1","theta=0.01","theta=0","allof=100","allof=50"),
       lty=c(1,1,1,1,1,2),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1),cex=0.6)

