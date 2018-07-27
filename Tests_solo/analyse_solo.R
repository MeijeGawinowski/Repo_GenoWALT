rm(list=ls())
setwd("~/Documents/GenoWALT/Tests_solo/Results")

ind=0:99
fitness_data <- data.frame(ind)
sum_vect <- c()

for (g in 0:20){
  file <- paste("pop",g,".csv",sep="")
  data <- read.table(file,sep="\t",header=TRUE)
  vect_fitness <- data$fitness
  sum_fitness <- sum(data$fitness)
  fitness_data <- cbind(fitness_data,vect_fitness)
  sum_vect <- c(sum_vect,sum_fitness)
}

for (g in 2:22){
  studied_fitness <- fitness_data[,g]
  hist(studied_fitness,main=paste("Histogramme des fitness à la génération",g-2,sep=""),
       ylab="Fréquence",xlab="Fitness",xlim=c(70,150),ylim=c(0,100))
}

plot(0:20,sum_vect,type="o",col="firebrick2",main="Evolution de la fitness globale",
     xlab="Génération",ylab="Fitness globale")
