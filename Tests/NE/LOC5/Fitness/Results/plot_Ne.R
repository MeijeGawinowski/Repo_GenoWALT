rm(list=ls())

setwd("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)

par(mfrow=c(1,1))

plot(dataA[,1],dataA[,2],col="indianred1",cex=0.8,ylim=c(0,1000))
for (i in seq(3:(dim(dataA)[2]-1))){
  points(dataA[,1],dataA[,i],col="indianred1",cex=0.8)  
}
abline(h=100,col="firebrick4")

plot(1:39,dataA$mean_Ne,col="indianred1",type="l")
abline(h=100,col="firebrick4",lty=2)

plot(dataB[,1],dataB[,2],col="lightgreen",cex=0.8,ylim=c(0,10000))
for (i in seq(3:(dim(dataB)[2]-1))){
  points(dataB[,1],dataB[,i],col="lightgreen",cex=0.8)  
}
abline(h=1000,col="forestgreen")

plot(1:99,dataB$mean_Ne,col="lightgreen",type="l")
abline(h=1000,col="forestgreen",lty=2)

plot(dataC[,1],dataC[,2],col="plum2",cex=0.8,ylim=c(0,500000))
for (i in seq(3:(dim(dataC)[2]-1))){
  print(i)
  points(dataC[,1],dataC[,i],col="plum2",cex=0.8)  
}
abline(h=10000,col="purple4")

plot(1:499,dataC$mean_Ne,col="plum2",type="l")
abline(h=10000,col="purple4",lty=2)
