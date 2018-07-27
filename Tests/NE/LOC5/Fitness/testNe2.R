rm(list=ls())
N_gen <- 40
time <- seq(1,N_gen-1)
tot_loc1 <- data.frame(time)
tot_loc2 <- data.frame(time)
tot_mloc <- data.frame(time)
N_rep <- 5

for (rep in seq(0,N_rep-1)){
  print(rep)
  path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/sizeC/Rep",rep,sep="")
  setwd(path_rep)
  name_file <- paste("datNe",rep,".csv",sep="")
  df <- read.table(name_file,header=TRUE)
  fc_loc1 <- as.numeric(df[,1])
  tot_loc1 <- cbind(tot_loc1,fc_loc1)
  fc_loc2 <- as.numeric(df[,2])
  tot_loc2 <- cbind(tot_loc2,fc_loc2)
  fc_mloc <- as.numeric(df[,3])
  tot_mloc <- cbind(tot_mloc,fc_mloc)
}

m1 <- c()
m2 <- c()
m3 <- c()
for(i in time){
  print(i)
  row1 <- as.numeric(tot_loc1[i,]) # Ne pour locus 1 au temps i pour tous les rep
  print(mean(row1))
  row2 <- as.numeric(tot_loc2[i,]) # Ne pour locus 2 au temps i pour tous les rep
  row3 <- as.numeric(tot_mloc[i,]) # Ne moyenné au temps i pour tous les rep
  m1 <- c(m1,mean(row1))
  m2 <- c(m2,mean(row2))
  m3 <- c(m3,mean(row3))
}

data_final <- data.frame(m1,m2,m3)
plot(1:39,data_final$m1)
