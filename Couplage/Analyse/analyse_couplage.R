### analyse population_data (fitness)


# change repertory accordy to studied repetition
rm(list=ls())
vect_fitness=c()

for (g in 0:19){
  file <- paste("pop",g,".csv",sep="")
  data <- read.table(file,sep=",",header=TRUE)
  hist(data$fitness,main=paste("Histogramme des fitness à la génération ",g,sep=""),ylab="Fréquence",xlab="Fitness")
  tot_fitness <- sum(data$fitness)
  vect_fitness <- c(vect_fitness,tot_fitness)
}

plot(0:19,vect_fitness,type="o",ylab="Fitness globale",xlab="Génération")



### analyse sim_scheme (phénotypes)

rm(list=ls())

## change directory to change population (i.e generation)
## LOOP

data=read.table("variable_geno_params.csv",header=TRUE,sep="\t")
names_genotypes <- strsplit(as.character(data$genotypes),",")[[1]]
count_genotypes <- strsplit(as.character(data$genotype_numbers),",")[[1]]
nb_geno <- length(names_genotypes)
population_PH <- c()
for (i in 1:nb_geno){
  geno <- names_genotypes[i]
  count <- count_genotypes[i]
  col_name <- paste("Param_PlHeight_",geno,sep="")
  PH_val <- data[col_name][[1]]
  PH_vect <- rep(PH_val,count)
  population_PH <- c(population_PH,PH_vect)
}
hist(population_PH)