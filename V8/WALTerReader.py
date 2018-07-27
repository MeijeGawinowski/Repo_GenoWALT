#
# WALTerReader.py $
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

import random as rd
import os
import numpy as np
import pandas as pd


class WALTerReader():
    
    def __init__(self):
        self.walter_file = "population_data.csv" 
        self.geno_names = "geno_conversion.csv" 
    
    
    def WALTer2GenoWALT(self):    
        file = open(self.walter_file,"rb")
        df=pd.read_csv(file,delimiter="\t",index_col=False)
        file.close()
        header = list(df)[3:np.shape(df)[1]]
        walter=np.array(df)
        print("header=",header)
        
        file = open(self.geno_names,"rb")
        df=pd.read_csv(file,delimiter="\t",index_col=False,header=None)
        file.close()
        geno=np.array(df)
        print("geno=",geno)
            
        N_geno=np.shape(geno)[0]
        N_ind=np.shape(walter)[0]
        
        fitness=list(walter[:,1])
        
        names_geno=list(geno[:,0])
        #print("names geno=",names_geno)
        dict = {}
        for ind in range(N_ind):
            name = walter[ind,2]
            geno_idx = names_geno.index(name)
            geno_dict = eval(geno[geno_idx,1])
            dict[ind] = geno_dict
        #print("DICT",dict)
        return([dict,header,fitness,walter])

        
        
    def setConversion(self,data): 
    # data is the QTL table
        N_qtl = np.shape(data)[0]
        tab = {}
        for qtl in range(N_qtl): # for each qtl in data_qtl
            values_qtl = data[qtl,3]
            # for each type qtl
            subdict = {}
            subdict["[0,0]"] = []
            subdict["[0,1]"] = []
            subdict["[1,0]"] = []
            subdict["[1,1]"] = []
            for k in range(len(values_qtl)):
                val_qtl = values_qtl[k] # we take each value of the qtl (i.e for each trait)
                subdict["[0,0]"].append(0)
                if data[qtl,5][k] == "complete":
                    subdict["[0,1]"].append(val_qtl)
                    subdict["[1,0]"].append(val_qtl)
                if data[qtl,5][k] == "partial":
                    subdict["[0,1]"].append(0.5*val_qtl)
                    subdict["[1,0]"].append(0.5*val_qtl)
                if data[qtl,5][k] == "recessive":
                    subdict["[0,1]"].append(0)
                    subdict["[1,0]"].append(0)
            
                subdict["[1,1]"].append(val_qtl)
            tab[qtl] = subdict
        return(tab)
        
    
    # def getGenoPop(self,data):
    # # data is the QTL table
    #     pheno = self.readFilePheno()
    #     N_qtl = np.shape(data)[0]
    #     N_ind = pheno[2]
    #     print("CHECK",N_qtl)
    #     try :
    #         file = open(self.fileGeno,"rb")
    #         listGeno = eval(file.read())
    #         file.close()
    #         if len(list(listGeno.values())[0]) != N_qtl:
    #             print("Warning : incorrect file, genotype randomly generated (wrong number of qtls)")
    #             listGeno = {}
    #             for ind in range(N_ind):
    #                 subdict = {}
    #                 for qtl in range(N_qtl):
    #                     subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
    #                 listGeno[ind] = subdict
    #     except FileNotFoundError :
    #         listGeno = {}
    #         for ind in range(N_ind):
    #             subdict = {}
    #             for qtl in range(N_qtl):
    #                 subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
    #             listGeno[ind] = subdict
    #         print("Warning :incorrect file, genotype randomly generated (wrong file or path name)")
    #     return(listGeno)
    
    def getAllelTab(self,listGeno):
    # listGeno is the dictionary list of the population genotype
        N_ind = len(listGeno)
        N_allel = len(listGeno[0])
        tab_out = np.zeros((N_ind,N_allel),dtype=object)
        for i in range(N_ind):
            tab_out[i,:] = list(listGeno[i].values())
        return(tab_out)
    
    def getQTLtab(self,popGeno,data):
        #popGeno is the population genotype dictionnary
        N_ind = len(popGeno)
        conv = self.setConversion(data)
        N_qtl = len(conv)
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_ind):
            qtl_val = []
            for j in range(N_qtl):
                allel = popGeno[i][j] # value of allele for individual i at qtl j
                al1=str(allel).replace(" ","")
                # conv[j] is the list of alleles and values for qtl j
                list_allel = list(conv[j].keys()) # list of alleles for qtl j
                list_valqtl = list(conv[j].values()) # list of qtl j values for each allele
                # al2 = "["+str(allel[1])+","+str(allel[4])+"]"
                al3 = "["+str(allel[0])+","+str(allel[1])+"]"
                try:
                    idx = list_allel.index(al1) # allele indice in the list
                except ValueError:
                    idx = list_allel.index(str(list(allel)))
                val = list_valqtl[idx] # qtl value matching this allele
                qtl_val.append(val)
            tab_out[i,:] = qtl_val
        return(tab_out)
    
    def dataWALTer(self,data):
        dict_genoPop = self.WALTer2GenoWALT()[0]
        header=self.WALTer2GenoWALT()[1]
        fitness=self.WALTer2GenoWALT()[2]
        walter=self.WALTer2GenoWALT()[3]
        genoPop = self.getAllelTab(dict_genoPop)
        conv = self.setConversion(data)
        # pheno = self.readFilePheno()
        tab = self.getQTLtab(genoPop,data)
        return([genoPop,conv,tab,header,fitness,walter])

