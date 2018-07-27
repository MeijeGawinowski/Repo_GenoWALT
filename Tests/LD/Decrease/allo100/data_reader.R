setwd("~/Documents/LD_Tests/Decrease/Control")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

plot(0:24,dataA$mean_r,type="o",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
points(0:24,dataB$mean_r,type="o",col="mediumorchid")
points(0:24,dataC$mean_r,type="o",col="dodgerblue4")
points(0:24,dataD$mean_r,type="o",col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0"),lty=c(1,1,1,1),
        col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen"))
