setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

plot(0:24,dataA$mean_r,type="o",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
points(0:24,dataB$mean_r,type="o",col="mediumorchid")
points(0:24,dataC$mean_r,type="o",col="dodgerblue4")
points(0:24,dataD$mean_r,type="o",col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0","théorique","observé"),lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,19),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1))

library(Metrics)

###########################
### Theoric DL decrease ###
###########################

DL_decrease <- function(r,x){
  return( (1-(2*r)/3)^x )
}

vect_r <- c(0.5,0.1,0.01,0)
xplot <- seq(0,24,0.1)
xtest <- seq(0,24)
data_plot <- data.frame(xplot)
data_test <- data.frame(xtest)

for (i in vect_r){
  yplot <- DL_decrease(i,xplot)
  ytest <- DL_decrease(i,xtest)
  data_plot <- cbind(data_plot,yplot)
  data_test <- cbind(data_test,ytest)
}

colnames(data_plot) <- colnames(data_test) <- c("gen","rA","rB","rC","rD")


###############
### Fitness ###
###############

# PLOT 1

plot(data_plot$gen,data_plot$rA,type="l",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
lines(data_plot$gen,data_plot$rB,type="l",col="mediumorchid")
lines(data_plot$gen,data_plot$rC,type="l",col="dodgerblue4")
lines(data_plot$gen,data_plot$rD,type="l",col="forestgreen")


setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

points(0:24,dataA$mean_r,col="firebrick2")
points(0:24,dataB$mean_r,col="mediumorchid")
points(0:24,dataC$mean_r,col="dodgerblue4")
points(0:24,dataD$mean_r,col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0","théorique","observé"),lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,19),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1))

# PLOT 2
plot(data_plot$gen,data_plot$rA,type="l",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
lines(data_plot$gen,data_plot$rB,type="l",col="mediumorchid")
lines(data_plot$gen,data_plot$rC,type="l",col="dodgerblue4")
lines(data_plot$gen,data_plot$rD,type="l",col="forestgreen")

points(0:24,dataA$vect_r1,col="firebrick2")
points(0:24,dataA$vect_r2,col="firebrick2")
points(0:24,dataA$vect_r3,col="firebrick2")
points(0:24,dataA$vect_r4,col="firebrick2")
points(0:24,dataA$vect_r5,col="firebrick2")

points(0:24,dataB$vect_r1,col="mediumorchid")
points(0:24,dataB$vect_r2,col="mediumorchid")
points(0:24,dataB$vect_r3,col="mediumorchid")
points(0:24,dataB$vect_r4,col="mediumorchid")
points(0:24,dataB$vect_r5,col="mediumorchid")

points(0:24,dataC$vect_r1,col="dodgerblue4")
points(0:24,dataC$vect_r2,col="dodgerblue4")
points(0:24,dataC$vect_r3,col="dodgerblue4")
points(0:24,dataC$vect_r4,col="dodgerblue4")
points(0:24,dataC$vect_r5,col="dodgerblue4")

points(0:24,dataD$vect_r1,col="forestgreen")
points(0:24,dataD$vect_r2,col="forestgreen")
points(0:24,dataD$vect_r3,col="forestgreen")
points(0:24,dataD$vect_r4,col="forestgreen")
points(0:24,dataD$vect_r5,col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0","théorique","observé"),lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,19),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1))

# Test 1 (pearson)
corA <- c(cor.test(data_test$rA,dataA$mean_r)$estimate,cor.test(data_test$rA,dataA$mean_r)$p.value)
corB <- c(cor.test(data_test$rB,dataB$mean_r)$estimate,cor.test(data_test$rB,dataB$mean_r)$p.value)
corC <- c(cor.test(data_test$rC,dataC$mean_r)$estimate,cor.test(data_test$rC,dataC$mean_r)$p.value)
corD <- c(cor.test(data_test$rD,dataD$mean_r)$estimate,cor.test(data_test$rD,dataD$mean_r)$p.value)

data_corF <- data.frame(corA,corB,corC,corD)
rownames(data_corF) <- c("r2","p-value")
setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50")
write.table(data_corF,"test_Fitness.csv")

# Test 2 (rmse)
rmseA <- rmse(actual=dataA$mean_r,predicted=data_test$rA)
rmseB <- rmse(actual=dataB$mean_r,predicted=data_test$rB)
rmseC <- rmse(actual=dataC$mean_r,predicted=data_test$rC)
rmseD <- rmse(actual=dataD$mean_r,predicted=data_test$rD)

data_rmseF <- data.frame(rmseA,rmseB,rmseC,rmseD)
colnames(data_rmseF) <- c("A","B","C","D")
setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50")
write.table(data_rmseF,"rmse_Fitness.csv")


###############
### Control ###
###############

# PLOT 1

plot(data_plot$gen,data_plot$rA,type="l",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
lines(data_plot$gen,data_plot$rB,type="l",col="mediumorchid")
lines(data_plot$gen,data_plot$rC,type="l",col="dodgerblue4")
lines(data_plot$gen,data_plot$rD,type="l",col="forestgreen")


setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Control")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

points(0:24,dataA$mean_r,col="firebrick2")
points(0:24,dataB$mean_r,col="mediumorchid")
points(0:24,dataC$mean_r,col="dodgerblue4")
points(0:24,dataD$mean_r,col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0","théorique","observé"),lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,19),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1))



# PLOT 2
plot(data_plot$gen,data_plot$rA,type="l",col="firebrick2",ylim=c(0,1.5),
     ylab="DL",xlab="t")
lines(data_plot$gen,data_plot$rB,type="l",col="mediumorchid")
lines(data_plot$gen,data_plot$rC,type="l",col="dodgerblue4")
lines(data_plot$gen,data_plot$rD,type="l",col="forestgreen")

points(0:24,dataA$vect_r1,col="firebrick2")
points(0:24,dataA$vect_r2,col="firebrick2")
points(0:24,dataA$vect_r3,col="firebrick2")
points(0:24,dataA$vect_r4,col="firebrick2")
points(0:24,dataA$vect_r5,col="firebrick2")

points(0:24,dataB$vect_r1,col="mediumorchid")
points(0:24,dataB$vect_r2,col="mediumorchid")
points(0:24,dataB$vect_r3,col="mediumorchid")
points(0:24,dataB$vect_r4,col="mediumorchid")
points(0:24,dataB$vect_r5,col="mediumorchid")

points(0:24,dataC$vect_r1,col="dodgerblue4")
points(0:24,dataC$vect_r2,col="dodgerblue4")
points(0:24,dataC$vect_r3,col="dodgerblue4")
points(0:24,dataC$vect_r4,col="dodgerblue4")
points(0:24,dataC$vect_r5,col="dodgerblue4")

points(0:24,dataD$vect_r1,col="forestgreen")
points(0:24,dataD$vect_r2,col="forestgreen")
points(0:24,dataD$vect_r3,col="forestgreen")
points(0:24,dataD$vect_r4,col="forestgreen")
points(0:24,dataD$vect_r5,col="forestgreen")
legend("topright",c("r=0.5","r=0.1","r=0.01","r=0","théorique","observé"),lty=c(1,1,1,1,1,NA),
       pch=c(NA,NA,NA,NA,NA,19),col=c("firebrick2","mediumorchid","dodgerblue4","forestgreen",1,1))



# TEST 1 (pearson)
corA <- c(cor.test(data_test$rA,dataA$mean_r)$estimate,cor.test(data_test$rA,dataA$mean_r)$p.value)
corB <- c(cor.test(data_test$rB,dataB$mean_r)$estimate,cor.test(data_test$rB,dataB$mean_r)$p.value)
corC <- c(cor.test(data_test$rC,dataC$mean_r)$estimate,cor.test(data_test$rC,dataC$mean_r)$p.value)
corD <- c(cor.test(data_test$rD,dataD$mean_r)$estimate,cor.test(data_test$rD,dataD$mean_r)$p.value)

data_corC <- data.frame(corA,corB,corC,corD)
rownames(data_corC) <- c("r2","p-value")

setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50")
write.table(data_corC,"test_Control.csv")

# Test 2 (rmse)
rmseA <- rmse(actual=dataA$mean_r,predicted=data_test$rA)
rmseB <- rmse(actual=dataB$mean_r,predicted=data_test$rB)
rmseC <- rmse(actual=dataC$mean_r,predicted=data_test$rC)
rmseD <- rmse(actual=dataD$mean_r,predicted=data_test$rD)

data_rmseC <- data.frame(rmseA,rmseB,rmseC,rmseD)
colnames(data_rmseC) <- c("A","B","C","D")
setwd("/home/meije/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50")
write.table(data_rmseC,"rmse_Control.csv")
