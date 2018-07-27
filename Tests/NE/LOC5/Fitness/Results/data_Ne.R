#################
### Functions ###
#################

calc_Ne <- function(tab,eff,ls.comb){
  lab<-tab[,1:2]
  
  comb=NULL
  for (i in 1:nrow(ls.comb))
  {
    pop.ini<-ls.comb[i,"pop_ini"]
    pop.fin<-ls.comb[i,"pop_fin"]
    frq.ini<-tab[,which(colnames(tab)==ls.comb[i,"pop_ini"])]
    frq.fin<-tab[,which(colnames(tab)==ls.comb[i,"pop_fin"])]
    comb.temp<-cbind.data.frame(lab,pop.ini,pop.fin,frq.ini,frq.fin)
    comb<-rbind.data.frame(comb,comb.temp)
  }
  
  svg<-comb    
  comb$fc.all<-(comb[,"frq.ini"]-comb[,"frq.fin"])^2/((comb[,"frq.ini"]+comb[,"frq.fin"])/2-(comb[,"frq.ini"]*comb[,"frq.fin"]))
  comb<-comb[!is.na(comb$fc.all),]
  comb$mgfac<-paste(comb[,c("marker")],comb[,c("pop.ini")],comb[,c("pop.fin")])
  fc.loc<-tapply(comb[,c("fc.all")],comb$mgfac,mean)
  nb.all<-tapply(comb[,"fc.all"],comb$mgfac,length)
  mg1<-merge(fc.loc,nb.all,by=0)
  fc.loc<-cbind.data.frame(t(data.frame(strsplit(mg1[,1]," "))),mg1[,-1])                      
  rownames(fc.loc)<-seq(1:nrow(fc.loc))
  colnames(fc.loc)<-c("marker","pop.ini","pop.fin","fc.loc","nb.all")    
  fc.loc$w.fc.loc<-fc.loc$fc.loc*fc.loc$nb.all
  fc.loc$mgfac<-paste(fc.loc[,c("pop.ini")],fc.loc[,c("pop.fin")])
  fc.loc$mgfac<-as.factor(as.character(fc.loc$mgfac))  
  sum.w.fc.loc<-tapply(fc.loc$w.fc.loc,fc.loc$mgfac,sum)
  sum.nb.all<-tapply(fc.loc$nb.all,fc.loc$mgfac,sum)
  mg2<-merge(sum.w.fc.loc,sum.nb.all,by=0)    
  fc.mean<-cbind.data.frame(t(data.frame(strsplit(mg2[,1]," "))),mg2[,-1])    
  rownames(fc.mean)<-seq(1,nrow(fc.mean))
  colnames(fc.mean)<-c("pop.ini","pop.fin","sum.w.fc.loc","sum.nb.all")
  fc.mean$mean.fc<-fc.mean$sum.w.fc.loc/fc.mean$sum.nb.all
  
  mg.ini<-merge(fc.mean,eff,by.x="pop.ini",by.y="pop_farm")
  colnames(mg.ini)[(ncol(mg.ini)-1):ncol(mg.ini)]<-c("generation.ini","demography.ini")
  mg.fin<-merge(mg.ini,eff,by.x="pop.fin",by.y="pop_farm")
  mg.fin$nb.generation<-mg.fin$generation-mg.fin$generation.ini
  mg.fin$Ne<-mg.fin$nb.generation/(2*mg.fin$mean.fc) #-1/(2*mg.fin$eff.y)-1/(2*mg.fin$eff.x))) 
  
  ddl<-mg.fin$sum.nb.all-17
  Xhi.inf<-qchisq(0.275,ddl)
  Xhi.sup<-qchisq(0.975,ddl)  
  mg.fin$ICr<-(ddl*mg.fin$mean.fc)/Xhi.inf
  mg.fin$ICl<-(ddl*mg.fin$mean.fc)/Xhi.sup
  
  mg.fin$Ne.ICl<-mg.fin$nb.generation/(2*(mg.fin$ICr-1/(2*mg.fin$eff.y)-1/(2*mg.fin$eff.x)))  
  mg.fin$Ne.ICr<-mg.fin$nb.generation/(2*(mg.fin$ICl-1/(2*mg.fin$eff.y)-1/(2*mg.fin$eff.x)))
  
  tab.NE.Fc<-cbind.data.frame(pop_ini=mg.fin$pop.ini,eff_ini=mg.fin$eff.x,demo_ini=mg.fin$demography.ini,
                              pop_fin=mg.fin$pop.fin,eff_fin=mg.fin$eff.y,demo_fin=mg.fin$demography,
                              nb_generation=mg.fin$nb.generation,Fc=mg.fin$mean.fc,Ne=mg.fin$Ne,ICleft=mg.fin$Ne.ICl,ICright=mg.fin$Ne.ICr)
  return(tab.NE.Fc)
}

rep_Ne <- function(N_gen, N_ind, N_rep,size){
  final_data <- data.frame(seq(1,N_gen-1))
  for (rep in seq(0,N_rep-1)){
    
    path_rep <- paste("~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/size",size,"/Rep",rep,sep="")
    setwd(path_rep)
  
    names <- c()
    for (gen in seq(0,N_gen-1)){
      name_pop <- paste("pop",gen,sep="")
      names <- c(names,name_pop)
    }
    
    v_eff <- demography <- rep(100,N_gen)
    generation <- seq(0,N_gen-1)
    
    eff <- data.frame(v_eff,names,generation,demography)
    colnames(eff) <- c("eff","pop_farm","generation","demography")
    pop_ini <- rep(names[1],N_gen-1)
    pop_fin <- names[2:N_gen]
    ls.comb <- data.frame(pop_ini,pop_fin)
    
    
    dat <- read.table("pop0.txt",header=TRUE)
    # loc 1
    loc1 <- c(dat[,1],dat[,2])
    freq10 <- length(which(loc1==0))/(2*N_ind)
    freq11 <- length(which(loc1==1))/(2*N_ind)
    # loc 2
    loc2 <- c(dat[,3],dat[,4])
    freq20 <- length(which(loc2==0))/(2*N_ind)
    freq21 <- length(which(loc2==1))/(2*N_ind)
    
    marker <- c("qPh1_0","qPh1_1","qPh2_0","qPh2_1")
    allele <- c(0,1,0,1)
    pop0 <- c(freq10,freq11,freq20,freq21)
    tab <- data.frame(marker,allele,pop0)
    
    for (gen in seq(1,N_gen-1)){
      file <- paste("pop",gen,".txt",sep="")
      dat <- read.table(file,header=TRUE)
      # loc 1
      loc1 <- c(dat[,1],dat[,2])
      freq10 <- length(which(loc1==0))/(2*N_ind)
      freq11 <- length(which(loc1==1))/(2*N_ind)
      # loc 2
      loc2 <- c(dat[,3],dat[,4])
      freq20 <- length(which(loc2==0))/(2*N_ind)
      freq21 <- length(which(loc2==1))/(2*N_ind)
      newcol <- c(freq10,freq11,freq20,freq21)
      names_tab <- colnames(tab)
      tab <- cbind(tab,newcol)
      colnames(tab) <- c(names_tab,paste("pop",gen,sep=""))
    }
    
    tabNe <- calc_Ne(tab,eff,ls.comb)
    
    v_gen <- tabNe$nb_generation
    v_Ne <- tabNe$Ne
    sorted_Ne <- c()
    for (i in seq(1,N_gen-1)){
      idx <- which(v_gen==i)
      Ne <- v_Ne[i]
      sorted_Ne <- c(sorted_Ne,Ne)
    }
    final_data <- cbind(final_data,sorted_Ne)
  }
  
  mean_Ne <- c()
  for (k in seq(1:dim(final_data)[1])){
    row <- as.numeric(final_data[k,])
    m <- mean(row)
    mean_Ne <- c(mean_Ne,m)
  }
  
  final_data <- cbind(final_data,mean_Ne)
  return(final_data)
}


##############
### size A ###
##############


dataA <- rep_Ne(N_gen=40,N_ind=100,N_rep=100,size="A")
write.table(dataA,"~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/tabA.csv")



##############
### size B ###
##############


dataB <- rep_Ne(N_gen=100,N_ind=1000,N_rep=100,size="B")
write.table(dataB,"~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/tabB.csv")



##############
### size C ###
##############


dataC <- rep_Ne(N_gen=500,N_ind=10000,N_rep=100,size="C")
write.table(dataC,"~/Documents/GenoWALT/V4_tests/Tests/NE/LOC2/Fitness/Results/tabC.csv")


