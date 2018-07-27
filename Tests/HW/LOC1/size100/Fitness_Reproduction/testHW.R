#######################
####### allo = 0 ######
#######################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/size100/Fitness_Reproduction/allo0")

pval <- c()
SE <- c()
WC <- c()
RH <- c()

N_gen <- 10
for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
  RH <- c(RH,as.numeric(row2[16]))
}


gen <- 0:9

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC,RH=RH)
write.table(data,"HW_allo0.csv")

#########################
####### allo = 0.5 ######
#########################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/size100/Fitness_Reproduction/allo0.5")
pval <- c()
SE <- c()
WC <- c()
RH <- c()

N_gen <- 10
for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
  RH <- c(RH,as.numeric(row2[16]))
}


gen <- 0:9

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC,RH=RH)
write.table(data,"HW_allo0.5.csv")

#######################
####### allo = 1 ######
#######################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/size100/Fitness_Reproduction/allo1")
pval <- c()
SE <- c()
WC <- c()
RH <- c()

N_gen <- 10
for (i in 0:9){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
  RH <- c(RH,as.numeric(row2[16]))
}


gen <- 0:9

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC,RH=RH)
write.table(data,"HW_allo1.csv")

###################
####### PLOT ######
###################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/size100/Fitness_Reproduction")

data1 <- read.table("allo0/HW_allo0.csv")
data2 <- read.table("allo0.5/HW_allo0.5.csv")
data3 <- read.table("allo1/HW_allo1.csv")

pdf("HW_Fitness_100.pdf")
par(mfrow=c(2,2))
plot(data1$gen,data1$pval,type="o",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=100, Fitness",ylim=c(0,2))
points(data1$gen,data2$pval, type="o",col="dodgerblue4")
points(data1$gen,data3$pval,type="o",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$SE,type="o",col="firebrick2",xlab="Generation",ylab="SE",main="Size=100, Fitness",ylim=c(0,0.01))
points(data1$gen,data2$SE, type="o",col="dodgerblue4")
points(data1$gen,data3$SE,type="o",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$WC,type="o",col="firebrick2",xlab="Generation",ylab="W&C Fis",main="Size=100, Fitness",ylim=c(-1.5,1))
points(data1$gen,data2$WC, type="o",col="dodgerblue4")
points(data1$gen,data3$WC,type="o",col="forestgreen")
legend("bottomright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$RH,type="o",col="firebrick2",xlab="Generation",ylab="R&H Fis",main="Size=100, Fitness",ylim=c(-1.5,1))
points(data1$gen,data2$RH, type="o",col="dodgerblue4")
points(data1$gen,data3$RH,type="o",col="forestgreen")
legend("bottomright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

dev.off()
