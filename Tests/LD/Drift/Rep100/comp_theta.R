########################
### Theta Comparison ###
########################

rm(list=ls())

# theta=0.5
setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")

dataA1 <- read.table("tabA.csv",header=TRUE)
dataB1 <- read.table("tabB.csv",header=TRUE)
dataC1 <- read.table("tabC.csv",header=TRUE)
dataD1 <- read.table("tabD.csv",header=TRUE)

# theta=0.01
setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/thetaC/Control")

dataA2 <- read.table("tabA.csv",header=TRUE)
dataB2 <- read.table("tabB.csv",header=TRUE)
dataC2 <- read.table("tabC.csv",header=TRUE)
dataD2 <- read.table("tabD.csv",header=TRUE)

par(mfrow=c(1,1))
plot(dataA1$gen,dataA1$mean_r,type="l",col="firebrick4",ylim=c(0,1),xlab="t",ylab="DL")
lines(dataA2$ge,dataA2$mean_r,type="l",col="firebrick4",lty=2)
lines(dataB1$ge,dataB1$mean_r,type="l",col="forestgreen")
lines(dataB2$ge,dataB2$mean_r,type="l",col="forestgreen",lty=2)
lines(dataC1$ge,dataC1$mean_r,type="l",col="purple4")
lines(dataC2$ge,dataC2$mean_r,type="l",col="purple4",lty=2)
lines(dataD1$gen,dataD1$mean_r,type="l",col="navy")
lines(dataD2$gen,dataD2$mean_r,type="l",col="navy",lty=2)
legend("topright",c("size=100","size=500","size=1000","size=10 000","theta=0.5","theta=0.01"),
       col=c("firebrick2","forestgreen","purple4","navy",1,1),lty=c(1,1,1,1,1,2))
