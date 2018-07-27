#######################
####### allo = 0 ######
#######################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen/allo0")

pval1 <- c()
SE1 <- c()
WC1 <- c()


for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
}

gen <- 0:99

data <- data.frame(gen,pval1,SE1,WC1)
write.table(data,"HW_allo0.csv")

#########################
####### allo = 0.5 ######
#########################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen/allo0.5")

pval1 <- c()
SE1 <- c()
WC1 <- c()


for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
}

gen <- 0:99

data <- data.frame(gen,pval1,SE1,WC1)
write.table(data,"HW_allo0.5.csv")

#######################
####### allo = 1 ######
#######################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen/allo1")

pval1 <- c()
SE1 <- c()
WC1 <- c()


for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row1 <- strsplit(readLines(res_name)[31]," ")[[1]]
  pval1 <- c(pval1,as.numeric(row1[9]))
  SE1 <- c(SE1,as.numeric(row1[11]))
  WC1 <- c(WC1,as.numeric(row1[14]))
}

gen <- 0:99

data <- data.frame(gen,pval1,SE1,WC1)
write.table(data,"HW_allo1.csv")


###################
####### PLOT ######
###################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/size4/100gen")

data1 <- read.table("allo0/HW_allo0.csv")
data2 <- read.table("allo0.5/HW_allo0.5.csv")
data3 <- read.table("allo1/HW_allo1.csv")

pdf("Fitness_size4.pdf")

par(mfrow=c(2,2))

plot(data1$gen,data1$pval1,type='l',main="Size 100 000, Fitness",ylab="p-value",
     xlab="Generations",col="firebrick2",ylim=c(0,1.3))
points(data1$gen,data2$pval1,type="l",col="dodgerblue4")
points(data1$gen,data3$pval1,type="l",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


plot(data1$gen,data1$SE1,type='l',main="Size 100 000, Fitness",ylab="SE",
     xlab="Generations",col="firebrick2",ylim=c(0,0.12))
points(data1$gen,data2$SE1,type="l",col="dodgerblue4")
points(data1$gen,data3$SE1,type="l",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data1$gen,data1$WC1,type='l',main="Size 100 000, Fitness",ylab="W&C Fis",
     xlab="Generations",col="firebrick2",ylim=c(-0.15,1.3))
points(data1$gen,data2$WC1,type="l",col="dodgerblue4")
points(data1$gen,data3$WC1,type="l",col="forestgreen")
legend("topright",cex=0.5,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

dev.off()

