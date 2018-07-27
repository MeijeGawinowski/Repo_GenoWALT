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
        

    def Creation(self,tab_qtl,tab_loci,genoPop,N_ind):

        # Genetic information
        list_chrom = list(np.add(tab_qtl[:,1],-1))
        list_pos = list(tab_qtl[:,2])
        list_names = list(tab_qtl[:,4])
        print(list_pos)
        pop = sim.Population(size=N_ind, loci=[0]*self.nchrom, ploidy=self.ploidy, infoFields=['fitness'])
        pop.addLoci(chrom=list(np.add(tab_loci[:,1],-1)), pos=list(tab_loci[:,2]))
        pop.addLoci(chrom=list_chrom, pos=list_pos, lociNames=list_names)
        
        # Information fields initialization
        fitness = [1/N_ind]*N_ind
        pop.setIndInfo(fitness,'fitness')
        
        # Retrieval of QTL indices for traits of interest
        qtl_names = tab_qtl[:,4]
        qtl_indices = list(pop.lociByNames(qtl_names))
        
        # Genotype initialization
        for i in range(N_ind):
            ind = pop.individual(i) # individual i
            val = genoPop[i,:] # individual i alleles
            for qtl in range(len(val)):
                allel = val[qtl]
                locus = qtl_indices[qtl]
                ind.setAllele(int(allel[0]),locus,0)
                ind.setAllele(int(allel[1]),locus,1)
        return(pop)


    


