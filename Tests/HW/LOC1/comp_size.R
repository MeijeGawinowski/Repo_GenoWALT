### FITNESS
rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1")

data1a <- read.table("size1000/Fitness_Reproduction/allo0/HW_allo0.csv")
data2a <- read.table("size1000/Fitness_Reproduction/allo0.5/HW_allo0.5.csv")
data3a <- read.table("size1000/Fitness_Reproduction/allo1/HW_allo1.csv")


data1b <- read.table("size3/Fitness_Reproduction/allo0/HW_allo0.csv")
data2b <- read.table("size3/Fitness_Reproduction/allo0.5/HW_allo0.5.csv")
data3b <- read.table("size3/Fitness_Reproduction/allo1/HW_allo1.csv")

data1c <- read.table("size4/100gen/HW_allo0.csv")
data2c <- read.table("size4/100gen/HW_allo0.5.csv")
data3c <- read.table("size4/100gen/HW_allo1.csv")


### p-value comparison
par(mfrow=c(1,3))
plot(data1a$gen,data1a$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=1000, Fitness",ylim=c(0,2))
points(data1a$gen,data2a$pval, type="l",col="dodgerblue4")
points(data1a$gen,data3a$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=10 000, Fitness",ylim=c(0,2))
points(data1b$gen,data2b$pval, type="l",col="dodgerblue4")
points(data1b$gen,data3b$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1c$gen,data1c$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=100 000, Fitness",ylim=c(0,2))
points(data1c$gen,data2c$pval, type="l",col="dodgerblue4")
points(data1c$gen,data3c$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

### SE comparison
par(mfrow=c(1,3))
plot(data1a$gen,data1a$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=1000, Fitness",ylim=c(0,0.03))
points(data1a$gen,data2a$SE, type="l",col="dodgerblue4")
points(data1a$gen,data3a$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=10 000, Fitness",ylim=c(0,0.1))
points(data1b$gen,data2b$SE, type="l",col="dodgerblue4")
points(data1b$gen,data3b$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1c$gen,data1c$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=100 000, Fitness",ylim=c(0,0.2))
points(data1c$gen,data2c$SE, type="l",col="dodgerblue4")
points(data1c$gen,data3c$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

### WC comparison
par(mfrow=c(1,3))
plot(data1a$gen,data1a$WC,type="l",col="firebrick2",xlab="Generation",ylab="WC",main="Size=1000, Fitness",ylim=c(0,1.5))
points(data1a$gen,data2a$WC, type="l",col="dodgerblue4")
points(data1a$gen,data3a$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$WC,type="l",col="firebrick2",xlab="Generation",ylab="WC",main="Size=10 000, Fitness",ylim=c(0,1.5))
points(data1b$gen,data2b$WC, type="l",col="dodgerblue4")
points(data1b$gen,data3b$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1c$gen,data1c$WC,type="l",col="firebrick2",xlab="Generation",ylab="WC",main="Size=100 000, Fitness",ylim=c(0,1.5))
points(data1c$gen,data2c$WC, type="l",col="dodgerblue4")
points(data1c$gen,data3c$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


### CONTROL
rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1")

data1a <- read.table("size1000/Controlled_Reproduction/allo0/HW_allo0.csv")
data2a <- read.table("size1000/Controlled_Reproduction/allo0.5/HW_allo0.5.csv")
data3a <- read.table("size1000/Controlled_Reproduction/allo1/HW_allo1.csv")


data1b <- read.table("size3/Controlled_Reproduction/allo0/HW_allo0.csv")
data2b <- read.table("size3/Controlled_Reproduction/allo0.5/HW_allo0.5.csv")
data3b <- read.table("size3/Controlled_Reproduction/allo1/HW_allo1.csv")


### p-value comparison
par(mfrow=c(1,2))
plot(data1a$gen,data1a$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=1000, Control",ylim=c(0,2))
points(data1a$gen,data2a$pval, type="l",col="dodgerblue4")
points(data1a$gen,data3a$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=10 000, Control",ylim=c(0,2))
points(data1b$gen,data2b$pval, type="l",col="dodgerblue4")
points(data1b$gen,data3b$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


### SE comparison
par(mfrow=c(1,2))
plot(data1a$gen,data1a$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=1000, Control",ylim=c(0,0.03))
points(data1a$gen,data2a$SE, type="l",col="dodgerblue4")
points(data1a$gen,data3a$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=10 000, Control",ylim=c(0,0.1))
points(data1b$gen,data2b$SE, type="l",col="dodgerblue4")
points(data1b$gen,data3b$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


### WC comparison
par(mfrow=c(1,2))
plot(data1a$gen,data1a$WC,type="l",col="firebrick2",xlab="Generation",ylab="WC",main="Size=1000, Control",ylim=c(0,1.5))
points(data1a$gen,data2a$WC, type="l",col="dodgerblue4")
points(data1a$gen,data3a$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1b$gen,data1b$WC,type="l",col="firebrick2",xlab="Generation",ylab="WC",main="Size=10 000, Control",ylim=c(0,1.5))
points(data1b$gen,data2b$WC, type="l",col="dodgerblue4")
points(data1b$gen,data3b$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

