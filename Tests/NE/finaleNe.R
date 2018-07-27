rm(list=ls())

Fc_loc <- function(freq10,freq1t,freq20,freq2t){
  fc1 <- ((freq10-freq1t)^2)/((freq10+freq1t)/2-freq10*freq1t)
  fc2 <- ((freq20-freq2t)^2)/((freq20+freq2t)/2-freq20*freq2t)
  Fc <- 0.5*(fc1+fc2)
  return(Fc)
}

calc_Ne1 <- function(Fc,t,S0,St){
  Ne <- t/(2*Fc)
  return(Ne)
}

calc_Ne2 <- function(Fc,t,S0,St){
  Ne <- (t-2)/(2*Fc)
  return(Ne)
}

calc_Ne3 <- function(Fc,t,S0,St){
  Ne <- t/(2*Fc-1/S0-1/St)
  return(Ne)
}

calc_Ne4 <- function(Fc,t,S0,St){
  Ne <- (t-2)/(2*Fc-1/S0-1/St)
  return(Ne)
}

test_Ne <- function(N_loc,N_rep,N_gen,N_ind,size,f_Ne){
  Ne_A <- c()
  m_fc <- c()
  for (rep in seq(0,N_rep-1)){
    print(rep)
    path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC",toString(N_loc),"/Fitness/Results/size",size,"/Rep",rep,sep="")
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
    Ne <- f_Ne(Fc=Mean_Fc,t=N_gen,S0=N_ind,St=N_ind)
    Ne_A <- c(Ne_A,Ne)
    m_fc <- c(m_fc,Mean_Fc)
  }
  l=list(v_Ne=Ne_A,m_Ne=mean(Ne_A),v_mFc=m_fc,mmFc=mean(m_fc))
  return(l)
}


##############
### 5 LOCI ###
##############

tab5_100 <- data.frame()
tab5_250 <- data.frame()
tab5_1000 <- data.frame()

# 1
testA <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne1)
res <- as.numeric(ci(testA$v_Ne,confidence=0.975))[1:3]
tab5_100 <- rbind(tab5_100,res)
testB <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne1)
res <- as.numeric(ci(testB$v_Ne,confidence=0.975))[1:3]
tab5_250 <- rbind(tab5_250,res)
testC <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne1)
res <- as.numeric(ci(testC$v_Ne,confidence=0.975))[1:3]
tab5_1000 <- rbind(tab5_1000,res)

# 2
testA <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne2)
res <- as.numeric(ci(testA$v_Ne,confidence=0.975))[1:3]
tab5_100 <- rbind(tab5_100,res)
testB <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne2)
res <- as.numeric(ci(testB$v_Ne,confidence=0.975))[1:3]
tab5_250 <- rbind(tab5_250,res)
testC <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne2)
res <- as.numeric(ci(testC$v_Ne,confidence=0.975))[1:3]
tab5_1000 <- rbind(tab5_1000,res)

# 3
testA <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne3)
res <- as.numeric(ci(testA$v_Ne,confidence=0.975))[1:3]
tab5_100 <- rbind(tab5_100,res)
testB <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne3)
res <- as.numeric(ci(testB$v_Ne,confidence=0.975))[1:3]
tab5_250 <- rbind(tab5_250,res)
testC <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne3)
res <- as.numeric(ci(testC$v_Ne,confidence=0.975))[1:3]
tab5_1000 <- rbind(tab5_1000,res)

# 4
testA <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=100,size="A",f_Ne=calc_Ne4)
res <- as.numeric(ci(testA$v_Ne,confidence=0.975))[1:3]
tab5_100 <- rbind(tab5_100,res)
testB <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=250,size="B",f_Ne=calc_Ne4)
res <- as.numeric(ci(testB$v_Ne,confidence=0.975))[1:3]
tab5_250 <- rbind(tab5_250,res)
testC <- test_Ne(N_loc=5,N_rep=1000,N_gen=10,N_ind=1000,size="C",f_Ne=calc_Ne4)
res <- as.numeric(ci(testC$v_Ne,confidence=0.975))[1:3]
tab5_1000 <- rbind(tab5_1000,res)





###############
### 10 LOCI ###
###############

testA1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=100,size="A")
testB1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=250,size="B")
testC1 <- test_Ne(N_loc=10,N_rep=1000,N_gen=10,N_ind=1000,size="C")



print(ci(testA1$v_Ne,confidence=0.975))
print(ci(testB1$v_Ne,confidence=0.975))
print(ci(testC1$v_Ne,confidence=0.975))


###############
### 30 LOCI ###
###############

testA2 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=100,size="A")
testB2 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=250,size="B")
testC2 <- test_Ne(N_loc=30,N_rep=1000,N_gen=10,N_ind=1000,size="C")


print(ci(testA2$v_Ne,confidence=0.975))
print(ci(testB2$v_Ne,confidence=0.975))
print(ci(testC2$v_Ne,confidence=0.975))

###############
### 50 LOCI ###
###############

testA3 <- test_Ne(N_loc=50,N_rep=1000,N_gen=10,N_ind=100,size="A")
testB3 <- test_Ne(N_loc=50,N_rep=1000,N_gen=10,N_ind=250,size="B")
testC3 <- test_Ne(N_loc=50,N_rep=1000,N_gen=10,N_ind=1000,size="C")


print(ci(testA3$v_Ne,confidence=0.975))
print(ci(testB3$v_Ne,confidence=0.975))
print(ci(testC3$v_Ne,confidence=0.975))


############
### PLOT ###
############


par(mfrow=c(3,2))
boxplot(testA$v_Ne,testA1$v_Ne,testA2$v_Ne,testA3$v_Ne,
        names=c("5 loci","10 loci","30 loci","50 loci"),
        main="Effectif efficace, N=100",col="indianred1")
boxplot(testA$v_Ne,testA1$v_Ne,testA2$v_Ne,ylim=c(0,300),
        names=c("5 loci","10 loci","30 loci"),main="Effectif efficace, N=100",
        col="indianred1")
abline(h=100,col="firebrick4")

boxplot(testB$v_Ne,testB1$v_Ne,testB2$v_Ne,testB3$v_Ne,
        names=c("5 loci","10 loci","30 loci","50 loci"),
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
