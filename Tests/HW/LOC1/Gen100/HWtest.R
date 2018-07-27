#######################
####### allo = 0 ######
#######################

rm(list=ls())
library(genepop)
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/Gen100/allo0")

pval <- c()
SE <- c()
WC <- c()

for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
}


gen <- 0:99

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC)
write.table(data,"HW_allo0.csv")

#########################
####### allo = 0.5 ######
#########################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/Gen100/allo0.5")
pval <- c()
SE <- c()
WC <- c()

for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
}


gen <- 0:99

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC)
write.table(data,"HW_allo0.5.csv")

#######################
####### allo = 1 ######
#######################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/Gen100/allo1")
pval <- c()
SE <- c()
WC <- c()


for (i in 0:99){
  name <- paste('pop',i,'.txt',sep="")
  print(name)
  res_name <- paste('pop',i,'.txt.D',sep="")
  test_HW(name,which='deficit',res_name)
  row <- readLines(res_name)[31]
  row2 <- strsplit(row," ")[[1]]
  pval <- c(pval,as.numeric(row2[9]))
  SE <- c(SE,as.numeric(row2[11]))
  WC <- c(WC,as.numeric(row2[14]))
}


gen <- 0:99

data <- data.frame(gen=gen,pval=pval,SE=SE,WC=WC)
write.table(data,"HW_allo1.csv")
# attention fichier non valide à corriger

###################
####### PLOT ######
###################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/HW/LOC1/Gen100")

data1 <- read.table("allo0/HW_allo0.csv")
data2 <- read.table("allo0.5/HW_allo0.5.csv")
data3 <- read.table("allo1/HW_allo1.csv")

pdf("HW_Fitness_1000.pdf")
par(mfrow=c(2,2))
plot(data2$gen,data1$pval,type="l",col="firebrick2",xlab="Generation",ylab="P-value",main="Size=1000, Fitness",ylim=c(0,2))
points(data2$gen,data2$pval, type="l",col="dodgerblue4")
points(data2$gen,data3$pval,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data2$gen,data1$SE,type="l",col="firebrick2",xlab="Generation",ylab="SE",main="Size=1000, Fitness",ylim=c(0,0.03))
points(data2$gen,data2$SE, type="l",col="dodgerblue4")
points(data2$gen,data3$SE,type="l",col="forestgreen")
legend("topright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))

plot(data2$gen,data1$WC,type="l",col="firebrick2",xlab="Generation",ylab="W&C Fis",main="Size=1000, Fitness",ylim=c(-1.5,1))
points(data2$gen,data2$WC, type="l",col="dodgerblue4")
points(data2$gen,data3$WC,type="l",col="forestgreen")
legend("bottomright",cex=0.75,c("allo=0","allo=0.5","allo=1"),col=c("firebrick2","dodgerblue4","forestgreen"),lty=c(1,1,1))


dev.off()
