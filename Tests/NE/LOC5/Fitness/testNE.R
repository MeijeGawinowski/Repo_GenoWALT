rm(list=ls())

# df_freq <- data.frame()

Fc_loc <- function(freq10,freq1t,freq20,freq2t){
  fc1 <- ((freq10-freq1t)^2)/((freq10+freq1t)/2-freq10*freq1t)
  fc2 <- ((freq20-freq2t)^2)/((freq20+freq2t)/2-freq20*freq2t)
  Fc <- 0.5*(fc1+fc2)
  return(Fc)
}

calc_Ne <- function(Fc,t){
  Ne <- t/(2*Fc)
  return(Ne)
}

calc_Ne2 <- function(Fc,t,S0,St){
  Ne <- t/(2*Fc-1/S0-1/St)
  return(Ne)
}

N_loc <- 5
N_rep <- 1000
N_gen <- 10
N_ind <- 100

Ne_A <- c()
m_fc <- c()

for (rep in seq(0,N_rep-1)){
  print(rep)
  path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/sizeA/Rep",rep,sep="")
  
  ## Calcul des fréquences alléliques à t=0
  setwd(path_rep)
  dat0 <- read.table("pop0.txt",header=TRUE)
  dat_freq0 <- c()
  for (i in seq(1,2*N_loc,2)){
    # pour chaque locus i
    vi00 <- dat0[,i]
    vi10 <- dat0[,i+1]
    vitot0 <- c(vi00,vi10)
    freq10 <- length(which(vitot0==0))/(2*N_ind) # fréquence allèle 1 à t=0
    freq20 <- length(which(vitot0==1))/(2*N_ind) # fréquence allèle 2 à t=0
    dat_freq0 <- c(dat_freq0,freq10,freq20)
    print(dat_freq0)
  }

  ## Calcul des fréquences alléliques à t = 10
  datt <- read.table("pop10.txt",header=TRUE)
  dat_freqt <- c()
  v_fc <- c()
  for (i in seq(1,2*N_loc,2)){
    # pour chaque locus i
    vi0t <- datt[,i]
    vi1t <- datt[,i+1]
    vitott <- c(vi0t,vi1t)
    freq1t <- length(which(vitott==0))/(2*N_ind) # fréquence allèle 1 à t=0
    freq2t <- length(which(vitott==1))/(2*N_ind) # fréquence allèle 2 à t=0
    dat_freqt <- c(dat_freqt,freq1t,freq2t)
    print(dat_freqt)
    fc <- Fc_loc(dat_freq0[i],freq1t,dat_freq0[i+1],freq2t)
    v_fc <- c(v_fc,fc)
    print(v_fc)
  }
  
  Mean_Fc <- mean(v_fc)
  Ne <- calc_Ne(Fc=Mean_Fc,t=N_gen-2)
  Ne_A <- c(Ne_A,Ne)
  m_fc <- c(m_fc,Mean_Fc)
}
  
  # for (gen in seq(1,N_gen-1)){
  #   v_fc <- c()
  #   # pour chaque génération > 0
  #   file <- paste("pop",gen,".txt",sep="")
  #   print(file)
  #   datt <- read.table(file,header=TRUE)
  #   dat_freqt <- c()
  #   for (i in seq(1,2*N_loc,2)){
  #     # pour chaque locus i
  #     print(paste("i",i,i+1))
  #     vi0t <- datt[,i]
  #     vi1t <- datt[,i+1]
  #     vitott <- c(vi0t,vi1t)
  #     freq1t <- length(which(vitott==0))/20000 # fréquence allèle 1 à t
  #     freq2t <- length(which(vitott==1))/20000 # fréquence allèle 2 à t
  #     print(freq1t,freq2t)
  #     dat_freqt <- c(dat_freqt,freq1t,freq2t)
  #     print(dat_freqt)
  #     fc <- Fc_loc(dat_freq0[i],freq1t,dat_freq0[i+1],freq2t)
  #     v_fc <- c(v_fc,fc)
  #   }
  #   print(v_fc)
  #   dat_fc <- rbind(dat_fc,v_fc)
  #   df_freq <- rbind(df_freq,dat_freqt)
  # }
  
  # mean_fc <- c()
  # for (k in 1:dim(dat_fc)[1]){
  #   row_dat <- as.numeric(dat_fc[k,])
  #   print(row_dat)
  #   m_fc <- mean(row_dat)
  #   mean_fc <- c(mean_fc,m_fc)
  # }
  # 
  # dat_fc <- cbind(dat_fc,mean_fc)
  # 
  # ## Calcul des Ne
  # time <- seq(1,dim(dat_fc)[1])
  # dat_Ne <- data.frame()
  # for (t in time){
  #   row_Fc <- as.numeric(dat_fc[t,]) # valeurs de Fc à chaque locus + moyenne au temps i
  #   print(row_Fc)
  #   row_Ne <- t/(2*row_Fc)
  #   print(row_Ne)
  #   dat_Ne <- rbind(dat_Ne,row_Ne)
  # } 
  # 
  # write.table(dat_fc,paste("datFc",rep,".csv",sep=""))
  # write.table(dat_Ne,paste("datNe",rep,".csv",sep=""))
  

# setwd("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/sizeB/Rep0")
# dat=read.table("datFc0.csv",header=TRUE)
# time=seq(1,39)
# mean_fc=dat$mean_fc/time
# datNe=read.table("datNe0.csv",header=TRUE)

print(Ne_A)
print(m_fc)
print(mean(Ne_A))
