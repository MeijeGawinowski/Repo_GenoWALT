rm(list=ls())
library(pegas)
setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/thetaD/Rep0")

vect_r1 <- c()

for (i in 0:24){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r1 <- abs(Mat[1])
  vect_r1 <- c(vect_r1,val_r1)
}

plot(0:24,vect_r1)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/thetaD/Rep1")
vect_r2 <- c()

for (i in 0:24){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r2 <- abs(Mat[1])
  vect_r2 <- c(vect_r2,val_r2)
}

plot(0:24,vect_r2)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/thetaD/Rep2")
vect_r3 <- c()

for (i in 0:24){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r3 <- abs(Mat[1])
  vect_r3 <- c(vect_r3,val_r3)
}

plot(0:24,vect_r3)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/thetaD/Rep3")
vect_r4 <- c()

for (i in 0:24){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r4 <- abs(Mat[1])
  vect_r4 <- c(vect_r4,val_r4)
}

plot(0:24,vect_r4)

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/thetaD/Rep4")
vect_r5 <- c()

for (i in 0:24){
  file=paste("pop",i,".txt",sep="")
  loc=read.loci(file)
  res=LD(loc)
  Mat=res$"Correlations among alleles"
  val_r5 <- abs(Mat[1])
  vect_r5 <- c(vect_r5,val_r5)
}

plot(0:24,vect_r5)

plot(0:24,vect_r1,type="l",col=1)
points(0:24,vect_r2,type="l",col=2)
points(0:24,vect_r3,type="l",col=3)
points(0:24,vect_r4,type="l",col=4)
points(0:24,vect_r5,type="l",col=5)

data_D <- data.frame(vect_r1,vect_r2,vect_r3,vect_r4,vect_r5)
mean_r <- c()
for (i in 1:dim(data_D)[1]){
  row_data <- as.numeric(data_D[i,])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}

newdata_D <- cbind(data_D,mean_r)
plot(0:24,newdata_D$mean_r)

write.table(newdata_D,"~/Documents/GenoWALT/V4_tests/Tests/LD/Decrease/allo50/Fitness/tabD.csv")
