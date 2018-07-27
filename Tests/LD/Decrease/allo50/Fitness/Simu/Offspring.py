#
# Offspring.py $/Users/meijegawinowski/Desktop/V4/Offspring.py
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
import csv
import copy as cp





class Offspring():


    def __init__(self,name):
        self.name = name


    def offspringGeno(self,offpop,conv,tab_qtl):
        """ This function returns the allel tab for the offspring population """
        # indices des QTL
        qtl_names = tab_qtl[:,4]
        qtl_indices = list(offpop.lociByNames(qtl_names))

        N_ind = offpop.popSize()
        N_qtl = len(qtl_indices)
        
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_qtl):
            list_allel = []
            for j in range(N_ind):
                ind = offpop.individual(j)
                val = [int(ind.allele(qtl_indices[i],0)),int(ind.allele(qtl_indices[i],1))] # for each individual j we take the QTL i and we take the allele value on the homoeologous copy 0 of the chromosome 
                list_allel.append(val)
            tab_out[:,i] = list_allel
        return(tab_out)
        
    def dictGenoOffspring(self,genotab):
        dict = {}
        for i in range(np.shape(genotab)[0]):
            subdict = {}
            for j in range(np.shape(genotab)[1]):
                subdict[j] = str(genotab[i,j])
            dict[i] = subdict
        return(dict)
       
        
        
    def Result(self,offpop,conv,tab_qtl):
        genotab = self.offspringGeno(offpop,conv,tab_qtl)
        geno_dict = self.dictGenoOffspring(genotab)
        return(genotab)
        

    

