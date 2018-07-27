#
# Population.py $
#
# This file is part of GenoWALT, a sofware of simulation of crop population.
#
# Copyright (C) 2017 Meije Gawinowski (meije.gawinowski@wanadoo.fr)
#
# 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#



import simuPOP as sim
import numpy as np


class Population():

    def __init__(self):
        self.nchrom = 21
        self.ploidy = 2
        self.trt_interest = 'PH'
        

    def Creation(self,tab_qtl,tab_loci,dict_qtl,genoPop,pheno):
        # Information on the population phenotype
        datPop = pheno[0] # phenotype table
        header = pheno[1] # traits
        N_ind = pheno[2] # population size

        # Genetic information
        list_chrom = list(np.add(tab_qtl[:,1],-1))
        list_pos = list(tab_qtl[:,2])
        list_names = list(tab_qtl[:,4])
        pop = sim.Population(size=N_ind, loci=[0]*self.nchrom, ploidy=self.ploidy, infoFields=header)
        pop.addLoci(chrom=list(np.add(tab_loci[:,1],-1)), pos=list(tab_loci[:,2]))
        pop.addLoci(chrom=list_chrom, pos=list_pos, lociNames=list_names)
        
        # Information fields initialization
        self.InfoPop(pop,header,datPop)
        
        # Retrieval of QTL indices
        interest_names = dict_qtl[self.trt_interest][3]
        idx_interest = list(pop.lociByNames(interest_names))
        
        # Genotype initialization
        for i in range(N_ind):
            ind = pop.individual(i) # individual i
            val = genoPop[i,:] # individual i alleles
            for k in range(len(val)):
                allel = val[k]
                locus = idx_interest[k]
                ind.setAllele(int(allel[0]),locus,0)
                ind.setAllele(int(allel[1]),locus,1)        
        return(pop)

    
    def InfoPop(self,pop,header,datPop):
        for i in range(len(header)):
            trt = header[i]
            val = list(datPop[:,i+1]) # vector of values for trait i
            pop.setIndInfo(val,trt)
            if trt == "offspring":
                offsp = val
        fitness = self.calcfitness(list(offsp))
        pop.addInfoFields('fitness')
        pop.setIndInfo(fitness,'fitness')
        return(pop)
        
    def calcfitness(self,offsp):
        all_offsp = sum(offsp)
        fitness = offsp/all_offsp
        return(list(fitness))

    def getMean(self,parpop):
        vect = parpop.indInfo(self.trt_interest)
        mean = np.mean(list(vect))
        return(mean)
    


