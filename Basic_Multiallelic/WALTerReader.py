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
import itertools as it


class WALTerReader():
    
    def __init__(self):
        self.fileGeno = "dict3.txt" # file with the population genotype
        self.filePheno = "data_pheno3.csv" # file with the population phenotype
        # self.dominance = "complete"
        self.distrib = "affine"
    
    def readFilePheno(self):
        file = open(self.filePheno,"rb")
        df = pd.read_csv(file,delimiter="\t")
        file.close()
        data = np.array(df)
        header = list(df)[1:np.shape(data)[1]]
        Nb = np.shape(data)[0]
        return([data,header,Nb])
    
    def affine(self,x,N,v):
        return((x*v)/(N-1))
    
    def expo(self,x,N,v):
        return((v*(math.exp(x)-1))/(math.exp(N-1)-1))
        
        
    def setConversion(self,data): 
        N_qtl = np.shape(data)[0]
        nb_allele = data[:,6]
        values_qtl = data[:,3] 
            
        tab={}
        for qtl in range(N_qtl):
            subdict = {}
            nb=nb_allele[qtl]
            
            if nb == 2:
                values_qtl = data[qtl,3]
                subdict["[0,0]"] = []
                subdict["[0,1]"] = []
                subdict["[1,0]"] = []
                subdict["[1,1]"] = []
                subdict["[0,0]"] = []
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
                    if data[qtl,5] == "recessive":
                        subdict["[0,1]"].append(0)
                        subdict["[1,0]"].append(0)
            
            else :
                # initialization of dictionary list for each arrangement
                values_list = data[qtl,3]
                arr=list(it.product(list(range(nb)),repeat=2))
                for k in range(len(arr)):
                    meanval_list = []
                    arr_elt = list(arr[k])
                    al1 = arr_elt[0]
                    al2 = arr_elt[1]
                    for j in range(len(values_list)):
                    # for v in values_list
                        v = values_list[j]
                        if self.distrib == "affine":
                            values=[self.affine(i,nb,v) for i in range(nb)] # allele values for qtl
                        if self.distrib == "expo":
                            values=[self.expo(i,nb,v) for i in range(nb)]
                        meanval = (values[al1]+values[al2])/2
                        meanval_list.append(meanval)
                    subdict[str(arr_elt)] = meanval
                
                    # for j in range(len(values)):
                    #     arr_elt = list(arr[k]) # list(tuple) avec les deux allèles
                    #     al1 = arr_elt[0]
                    #     al2 = arr_elt[1]
                    #     meanval = (values[al1]+values[al2])/2
                    #     subdict[str(arr_elt)].append(meanval)
            tab[qtl]=subdict
         
        return(tab)
        
    
    def getGenoPop(self,data):
    # data is the QTL table
        pheno = self.readFilePheno()
        N_qtl = np.shape(data)[0]
        N_ind = pheno[2]
        try :
            file = open(self.fileGeno,"rb")
            listGeno = eval(file.read())
            file.close()
            if len(list(listGeno.values())[0]) != N_qtl:
                print("Warning : incorrect file, genotype randomly generated (wrong number of qtls)")
                listGeno = {}
                for ind in range(N_ind):
                    subdict = {}
                    for qtl in range(N_qtl):
                        subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
                    listGeno[ind] = subdict
        except FileNotFoundError :
            listGeno = {}
            for ind in range(N_ind):
                subdict = {}
                for qtl in range(N_qtl):
                    subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
                listGeno[ind] = subdict
            print("Warning :incorrect file, genotype randomly generated (wrong file or path name)")
        return(listGeno)
    
    def getAllelTab(self,listGeno):
    # listGeno is the dictionary list of the population genotype
        N_ind = len(listGeno)
        N_allel = len(listGeno[0])
        tab_out = np.zeros((N_ind,N_allel),dtype=object)
        for i in range(N_ind):
            tab_out[i,:] = list(listGeno[i].values())
        return(tab_out)
    
    def getQTLtab(self,popGeno,data):
        N_ind = len(listGeno)
        conv = self.setConversion(data)
        N_qtl = len(conv)
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_ind):
            qtl_val = []
            for j in range(N_qtl):
                allel = listGeno[i][j] # value of allele i for qtl j
                # conv[j] is the list of alleles and values for qtl j
                list_allel = list(conv[j].keys()) # list of alleles for qtl j
                list_valqtl = list(conv[j].values()) # list of qtl j values for each allele
                idx = list_allel.index(allel) # allele indice in the list
                val = list_valqtl[idx] # qtl value matching this allele
                qtl_val.append(val)
            tab_out[i,:] = qtl_val
        return(tab_out)
    
    def dataPheno(self,data):
        dict_genoPop = self.getGenoPop(data)
        genoPop = self.getAllelTab(dict_genoPop)
        conv = self.setConversion(data)
        pheno = self.readFilePheno()
        return([genoPop,conv,pheno])

