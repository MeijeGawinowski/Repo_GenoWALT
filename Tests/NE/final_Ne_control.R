rm(list=ls())
library(gmodels)
Fc_loc <- function(freq10,freq1t,freq20,freq2t){
  fc1 <- ((freq10-freq1t)^2)/((freq10+freq1t)/2-freq10*freq1t)
  fc2 <- ((freq20-freq2t)^2)/((freq20+freq2t)/2-freq20*freq2t)
  Fc <- 0.5*(fc1+fc2)
  return(Fc)
}

calc_Ne1 <- function(Fc,t,N_ind){
  Ne <- t/(2*Fc)
  return(Ne)
}

calc_Ne2 <- function(Fc,t,N_ind){
  Ne <- (t-2)/(2*(Fc-1/N_ind))
  return(Ne)
}

test_Ne <- function(N_loc,N_rep,N_gen,N_ind,size,f_Ne){
  Ne_A <- c()
  m_fc <- c()
  for (rep in seq(0,N_rep-1)){
    print(rep)
    path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC",toString(N_loc),"/Control/Results/size",size,"/Rep",rep,sep="")
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
      fc <- Fc_loc(dat_freq0[i],freq1t,dat_freq0[i+1],freq2t)
      v_fc <- c(v_fc,fc)
    }
    
    Mean_Fc <- mean(v_fc)
    Ne <- f_Ne(Fc=Mean_Fc,t=N_gen,N_ind=N_ind)
    Ne_A <- c(Ne_A,Ne)
    m_fc <- c(m_fc,Mean_Fc)
  }
  l=list(v_Ne=Ne_A,m_Ne=mean(Ne_A),v_mFc=m_fc,mmFc=mean(m_fc))
  true_Ne <- calc_Ne1(mean(m_fc),N_gen,N_ind)
  return(true_Ne)
}


##############
### 5 LOCI ###
##############

testA_1 <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne1)

testB_1 <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne1)

testC_1 <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne1)


###############
### 10 LOCI ###
###############

testA1_1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne1)

testB1_1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne1)

testC1_1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne1)



###############
### 30 LOCI ###
###############

testA2_1 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne1)

testB2_1 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne1)

testC2_1 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne1)




############
### PLOT ###
############


par(mfrow=c(3,2))
boxplot(testA$v_Ne,testA1$v_Ne,testA2$v_Ne,names=c("5 loci","10 loci","30 loci"),
        main="Effectif efficace, N=100",col="indianred1")
boxplot(testA$v_Ne,testA1$v_Ne,testA2$v_Ne,ylim=c(0,300),
        names=c("5 loci","10 loci","30 loci"),main="Effectif efficace, N=100",
        col="indianred1")
abline(h=100,col="firebrick4")

boxplot(testB$v_Ne,testB1$v_Ne,testB2$v_Ne,names=c("5 loci","10 loci","30 loci"),
        main="Effectif efficace, N=250",col="lightgreen")
boxplot(testB$v_Ne,testB1$v_Ne,testB2$v_Ne,ylim=c(0,800),
        names=c("5 loci","10 loci","30 loci"),main="Effectif efficace, N=250",
        col="lightgreen")
abline(h=250,col="forestgreen")

boxplot(testC$v_Ne,testC1$v_Ne,testC2$v_Ne,names=c("5 loci","10 loci","30 loci"),
        main="Effectif efficace, N=1000",col="plum2")
boxplot(testC$v_Ne,testC1$v_Ne,testC2$v_Ne,ylim=c(0,4000),
        names=c("5 loci","10 loci","30 loci"),main="Effectif efficace, N=1000",
        col="plum2")
abline(h=1000,col="purple4")
