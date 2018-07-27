library(pegas)
####################
### size A (100) ###
####################

rm(list=ls())

main_path="~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control/sizeA/Results/"

gen <- 0:99
dataA <- data.frame(gen)




for (rep in 0:99){
  # pour chaque réplicat rep
  name_rep <- paste("Rep",rep,sep="")
  path_rep <- paste(main_path,"Rep",rep,sep="")
  setwd(path_rep) # on se place dans le répertoire rep
  
  vect <- c()
  
  for (i in 0:99){
    file=paste("pop",i,".txt",sep="")
    loc=read.loci(file)
    res=LD(loc)
    Mat=res$"Correlations among alleles"
    val <- abs(Mat[1])
    vect <- c(vect,val)
  }
  
  names_data <- colnames(dataA) # stock colnames
  dataA <- cbind(dataA,vect)
  colnames(dataA) <- c(names_data,paste("vect",rep,sep=""))
  
} 

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")
write.table(dataA,"tabA.csv")


####################
### size B (500) ###
####################

rm(list=ls())

main_path="~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control/sizeB/Results/"

gen <- 0:99
dataB <- data.frame(gen)



for (rep in 0:99){
  # pour chaque réplicat rep
  name_rep <- paste("Rep",rep,sep="")
  path_rep <- paste(main_path,"Rep",rep,sep="")
  setwd(path_rep) # on se place dans le répertoire rep
  
  vect <- c()
  
  for (i in 0:99){
    file=paste("pop",i,".txt",sep="")
    loc=read.loci(file)
    res=LD(loc)
    Mat=res$"Correlations among alleles"
    val <- abs(Mat[1])
    vect <- c(vect,val)
  }
  
  names_data <- colnames(dataB) # stock colnames
  dataB <- cbind(dataB,vect)
  colnames(dataB) <- c(names_data,paste("vect",rep,sep=""))
  
} 

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")
write.table(dataB,"tabB.csv")

#####################
### size C (1000) ###
#####################

rm(list=ls())

main_path="~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control/sizeC/Results/"

gen <- 0:99
dataC <- data.frame(gen)



for (rep in 0:99){
  # pour chaque réplicat rep
  name_rep <- paste("Rep",rep,sep="")
  path_rep <- paste(main_path,"Rep",rep,sep="")
  setwd(path_rep) # on se place dans le répertoire rep
  
  vect <- c()
  
  for (i in 0:99){
    file=paste("pop",i,".txt",sep="")
    loc=read.loci(file)
    res=LD(loc)
    Mat=res$"Correlations among alleles"
    val <- abs(Mat[1])
    vect <- c(vect,val)
  }
  
  names_data <- colnames(dataC) # stock colnames
  dataC <- cbind(dataC,vect)
  colnames(dataC) <- c(names_data,paste("vect",rep,sep=""))
  
} 

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")
write.table(dataC,"tabC.csv")

#######################
### size D (10 000) ###
#######################

rm(list=ls())

main_path="~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control/sizeD/Results/"

gen <- 0:99
dataD <- data.frame(gen)



for (rep in 0:99){
  # pour chaque réplicat rep
  name_rep <- paste("Rep",rep,sep="")
  path_rep <- paste(main_path,"Rep",rep,sep="")
  setwd(path_rep) # on se place dans le répertoire rep
  
  vect <- c()
  
  for (i in 0:99){
    file=paste("pop",i,".txt",sep="")
    loc=read.loci(file)
    res=LD(loc)
    Mat=res$"Correlations among alleles"
    val <- abs(Mat[1])
    vect <- c(vect,val)
  }
  
  names_data <- colnames(dataD) # stock colnames
  dataD <- cbind(dataD,vect)
  colnames(dataD) <- c(names_data,paste("vect",rep,sep=""))
  
} 

setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")
write.table(dataD,"tabD.csv")


##################
### Comparison ###
##################

rm(list=ls())
setwd("~/Documents/GenoWALT/V4_tests/Tests/LD/Drift/Rep100/Control")

dataA <- read.table("tabA.csv",header=TRUE)
dataB <- read.table("tabB.csv",header=TRUE)
dataC <- read.table("tabC.csv",header=TRUE)
dataD <- read.table("tabD.csv",header=TRUE)

## add mean column

# dataA
dataA2 <- dataA
dataA2 <- na.omit(dataA2)
mean_r <- c()
for (i in 1:dim(dataA)[1]){
  row_data <- as.numeric(dataA[i,2:dim(dataA)[2]])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}
dataA <- cbind(dataA,mean_r)

mean_r <- c()
for (i in 1:dim(dataA2)[1]){
  row_data <- as.numeric(dataA2[i,2:dim(dataA2)[2]])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}
dataA2 <- cbind(dataA2,mean_r)

mean_r <- c()
for (i in 1:dim(dataB)[1]){
  row_data <- as.numeric(dataB[i,2:dim(dataB)[2]])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}
dataB <- cbind(dataB,mean_r)

mean_r <- c()
for (i in 1:dim(dataC)[1]){
  row_data <- as.numeric(dataC[i,2:dim(dataC)[2]])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}
dataC <- cbind(dataC,mean_r)

mean_r <- c()
for (i in 1:dim(dataD)[1]){
  row_data <- as.numeric(dataD[i,2:dim(dataD)[2]])
  m_r <- mean(row_data)
  mean_r <- c(mean_r,m_r)
}
dataD <- cbind(dataD,mean_r)


### Plot
par(mfrow=c(2,3))
# dataA
plot(dataA$gen,dataA[,2],col="indianred1",ylim=c(0,1),cex=0.8,ylab="DL",xlab="t",
     main="Control, size=100, non corrigé")
for (i in (3:(dim(dataA)[2]-1))){
  print(i)
  y <- dataA[,i]
  points(dataA$gen,y,col="indianred1",cex=0.8)
}
lines(dataA$gen,dataA$mean_r,type="l",col="firebrick4")

plot(dataA2$gen,dataA2[,2],col="indianred1",ylim=c(0,1),cex=0.8,ylab="DL",xlab="t",main="Control, size=100, corrigé")
for (i in (3:(dim(dataA2)[2]-1))){
  print(i)
  y <- dataA2[,i]
  points(dataA2$gen,y,col="indianred1",cex=0.8)
}
lines(dataA2$gen,dataA2$mean_r,type="l",col="firebrick4")

# dataB
plot(dataB$gen,dataB[,2],col="lightgreen",ylim=c(0,1),cex=0.8,ylab="DL",xlab="t",main="Control, size=500")
for (i in (2:(dim(dataB)[2]-1))){
  print(i)
  y <- dataB[,i]
  points(dataB$gen,y,col="lightgreen",cex=0.8)
}
lines(dataB$gen,dataB$mean_r,type="l",col="forestgreen")

# dataC
plot(dataC$gen,dataC[,2],col="plum2",ylim=c(0,1),cex=0.8,ylab="DL",xlab="t",main="Control, size=1000")
for (i in (2:(dim(dataC)[2]-1))){
  print(i)
  y <- dataC[,i]
  points(dataC$gen,y,col="plum2",cex=0.8)
}
lines(dataC$gen,dataC$mean_r,type="l",col="purple4")

# dataD
plot(dataD$gen,dataD[,2],col="skyblue",ylim=c(0,1),cex=0.8,ylab="DL",xlab="t",main="Control, size=10 000")
for (i in (2:(dim(dataD)[2]-1))){
  print(i)
  y <- dataD[,i]
  points(dataD$gen,y,col="skyblue",cex=0.8)
}
lines(dataD$gen,dataD$mean_r,col="navy",type="l")

plot(dataA$gen,dataA$mean_r,col="firebrick2",ylim=c(0,0.6),type="l",ylab="DL",xlab="t",main="Control")
lines(dataB$gen,dataB$mean_r,col="forestgreen",type="l")
lines(dataC$gen,dataC$mean_r,col="purple4",type="l")
lines(dataD$gen,dataD$mean_r,col="navy",type="l")
legend("topright",c("size=100","size=500","size=1000","size=10 000"),cex=0.7,col=c("firebrick2","forestgreen","purple4","navy"),lty=c(1,1,1,1))
