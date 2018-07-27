# pvalue plot comparison
par(mfrow=c(1,3))
# 1000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size1000/Fitness_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$pval,type="o",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=1000, Fitness",ylim=c(0,2))
points(data11$gen,data21$pval, type="o",col="dodgerblue4")
points(data11$gen,data31$pval,type="o",col="forestgreen")
abline(h=0.05,lty=2)
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 10 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size3/Fitness_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=10 000, Fitness",ylim=c(0,2))
points(data11$gen,data21$pval, type="l",col="dodgerblue4")
points(data11$gen,data31$pval,type="l",col="forestgreen")
abline(h=0.05,lty=2)
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 100 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen")

data11 <- read.table("HW_allo0.csv")
data21 <- read.table("HW_allo0.5.csv")
data31 <- read.table("HW_allo1.csv")

plot(data11$gen,data11$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=100 000, Fitness",ylim=c(0,2))
points(data11$gen,data21$pval, type="l",col="dodgerblue4")
points(data11$gen,data31$pval,type="l",col="forestgreen")
abline(h=0.05,lty=2)
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


# pvalue plot comparison
par(mfrow=c(1,2))
# 1000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size1000/Controlled_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$pval,type="o",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=1000, Control",ylim=c(0,1))
points(data11$gen,data21$pval, type="o",col="dodgerblue4")
points(data11$gen,data31$pval,type="o",col="forestgreen")
abline(h=0.05,lty=2)
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 10 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size3/Controlled_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=10 000, Control",ylim=c(0,2))
points(data11$gen,data21$pval, type="l",col="dodgerblue4")
points(data11$gen,data31$pval,type="l",col="forestgreen")
abline(h=0.05,lty=2)
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))





# WC plot comparison
par(mfrow=c(1,3))
# 1000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size1000/Fitness_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$WC,type="o",col="firebrick2",xlab="Generation",ylab="FIS",main="Size=1000, Fitness",ylim=c(0,2))
points(data11$gen,data21$WC, type="o",col="dodgerblue4")
points(data11$gen,data31$WC,type="o",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 10 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size3/Fitness_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$WC,type="l",col="firebrick2",xlab="Generation",ylab="FIS",main="Size=10 000, Fitness",ylim=c(0,2))
points(data11$gen,data21$WC, type="l",col="dodgerblue4")
points(data11$gen,data31$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 100 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen")

data11 <- read.table("HW_allo0.csv")
data21 <- read.table("HW_allo0.5.csv")
data31 <- read.table("HW_allo1.csv")

plot(data11$gen,data11$WC,type="l",col="firebrick2",xlab="Generation",ylab="FIS",main="Size=100 000, Fitness",ylim=c(0,2))
points(data11$gen,data21$WC, type="l",col="dodgerblue4")
points(data11$gen,data31$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


# WC plot comparison
par(mfrow=c(1,2))
# 1000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size1000/Controlled_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$WC,type="o",col="firebrick2",xlab="Generation",ylab="FIS",main="Size=1000, Control",ylim=c(0,1))
points(data11$gen,data21$WC, type="o",col="dodgerblue4")
points(data11$gen,data31$WC,type="o",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

# 10 000
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size3/Controlled_Reproduction")

data11 <- read.table("allo0/HW_allo0.csv")
data21 <- read.table("allo0.5/HW_allo0.5.csv")
data31 <- read.table("allo1/HW_allo1.csv")

plot(data11$gen,data11$WC,type="l",col="firebrick2",xlab="Generation",ylab="FIS",main="Size=10 000, Control",ylim=c(0,2))
points(data11$gen,data21$WC, type="l",col="dodgerblue4")
points(data11$gen,data31$WC,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))
