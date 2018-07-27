#
# Offspring.py $
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





class Offspring():


    def __init__(self):
        self.trt_interest = 'PH'


    def offspringGeno(self,offpop,conv,dict_qtl):
        """ This function returns the allel tab for the offspring population """
        # indices des QTL
        interest_names = dict_qtl[self.trt_interest][3]
        idx_interest = list(offpop.lociByNames(interest_names))

        N_ind = offpop.popSize()
        N_qtl = len(idx_interest)
        
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_qtl):
            list_allel = []
            for j in range(N_ind):
                ind = offpop.individual(j)
                val = [int(ind.allele(idx_interest[i],0)),int(ind.allele(idx_interest[i],1))] # for each individual j we take the QTL i and we take the allele value on the homoeologous copy 0 of the chromosome 
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

    def offspringQTLTab(self,dictGeno,conv):
        N_ind = len(dictGeno)
        N_qtl = len(conv)
        tab_out = np.zeros((N_ind,N_qtl))
        for i in range(N_ind):
            qtl_val = []
            for j in range(N_qtl):
                allel = dictGeno[i][j] # allel value for the individual i at the qtl j       
                al2="["+str(allel[1])+","+str(allel[4])+"]"
                # conv[j] est liste allèles/valeurs pour le qtl j
                list_allel = list(conv[j].keys()) # alleles list for qtl j
                list_valqtl = list(conv[j].values()) # qtl j values list for each allele
                try :
                	idx = list_allel.index(allel) # index of the allele in the list
                except ValueError :
                	idx = list_allel.index(al2)
                val = list_valqtl[idx] # qtl value matching this allele
                qtl_val.append(val)
            tab_out[i,:] = qtl_val
        return(tab_out)
    
    def rescaling(self,init,real_min,real_max):
        init_min = min(init)
        if init_min < 0:
            init = list(np.add(-init_min,init))
        y = [(real_max-real_min)/(max(init)-min(init))*init[i]+real_min for i in range(len(init))]
        return(y)
    
    
    
    def offspringPheno(self,qtltab,trt_mean,choice):
        N_ind = np.shape(qtltab)[0]
        tab_pheno = np.zeros((N_ind,2))
        if choice == "file" :
            PH = [trt_mean+sum(qtltab[i,:]) for i in range(N_ind)]
        if choice == "sampling" :
            PH = [sum(qtltab[i,:]) for i in range(N_ind)]
            PH = self.rescaling(PH,57.9,121.6)
        tab_pheno[:,0] = range(N_ind)
        tab_pheno[:,1] = PH
        return(tab_pheno)
        
        
        
    def Result(self,offpop,conv,dict_qtl,trt_mean,choice):
        genotab = self.offspringGeno(offpop,conv,dict_qtl)
        geno_dict = self.dictGenoOffspring(genotab)
        qtltab = self.offspringQTLTab(geno_dict,conv)
        tabpheno = self.offspringPheno(qtltab,trt_mean,choice)
        
        # Writing of a genotype dictionary in a txt file
        res_dict_file = open("dictRes.txt","w")
        res_dict_file.write(str(geno_dict))
        res_dict_file.close()
        
        # Writing of a PH phenotype table in a csv file
        res_pheno_file = open("phenoRes.csv","w")
        pheno_writer = csv.writer(res_pheno_file,delimiter="\t")
        pheno_writer.writerows(tabpheno)
        res_pheno_file.close()
       
        return([geno_dict,tabpheno])
        

    

