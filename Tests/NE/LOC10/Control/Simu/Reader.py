#
# Reader.py $
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
 

import pandas as pd
import numpy as np
import os



class Reader:
    
    def __init__(self):
        self.fileQTL = "dataQTL_D.csv" # file with genetic information of wheat
        self.choice = "file" # choice=file will take filteQTL to generate QTL data and choice="sampling" will generate loci and qtl randomly distributed in the genome
        self.met = "normal" # if choice=sampling and if met="normal" qtl values are sampled in a normal distribution and if met="exponential" they are sampled in an exponential distribution
        self.list_trt = ['PH']

    
    def dataQTL(self):
        if self.choice == "file":
            print(self.fileQTL)
            file = open(self.fileQTL,"rb")
            df=pd.read_csv(file,delimiter=",",index_col=False)
            file.close()
            print("df",df)
            data = np.array(df)
            qtl_indices = []
            for trt in self.list_trt :
                qtl_indices = qtl_indices + self.all_indices(trt,list(data[:,0]))
            loci_indices = self.all_indices('locus',list(data[:,0]))
            data_qtl = data[qtl_indices,:]
            data_loci = data[loci_indices,:]
            data_loci = np.delete(data_loci,[3,4],1)
            
        if self.choice == "sampling":
            N_qtl = int(input("Nombre de qtl : "))
            N_loci = int(input("Nombre de loci : "))
            type_qtl = list(np.random.choice(self.list_trt,size=N_qtl,replace=True))
            
            N_tot = N_loci+N_qtl
        
            chrom_qtl = np.random.choice(range(1,22),N_qtl,replace=True)
            pos_qtl = list(np.random.uniform(30,200,N_qtl))
            if self.met == "normal":
                val = list(np.random.normal(loc=0,scale=1,size=N_qtl))
            if self.met == "exponential":
                val = np.random.exponential(scale=1,size=N_qtl)
            names = ['qPh'+str(i) for i in range(N_qtl)]
            dominance = np.random.choice(['complete','partial','recessive'],size=N_qtl,replace=True)
            data_qtl = np.column_stack([type_qtl,chrom_qtl.astype(np.object),pos_qtl,val,names,dominance])
            type_loci = ['locus']*N_loci
            chrom_loci = np.random.choice(range(1,22),N_loci,replace=True)
            pos_loci = list(np.random.uniform(30,200,N_loci))
            data_loci = np.column_stack([type_loci,chrom_loci.astype(np.object),pos_loci])
        
        qtl = self.pleiotropyManagement(data_qtl)
        return([qtl, data_loci])
    


    def all_indices(self,value,qlist):
        """ This function returns a list with the indexes of all the occurences of a given value in a given list. """
        indices = []
        idx = -1
        while True:
            try:
                idx = qlist.index(value,idx+1)
                indices.append(idx)
            except ValueError:
                break
        return(indices)
    
        
    def dictGenerator(self,data):
        """ This function generates a dictionnary list with each trait as key and list of chromosomes, list of loci positions and list of loci names as values. Pleiotropy isn't taken into account. """
        qtl_keys = list(np.unique(data[:,0])) # unique values of traits
        dic = {}
        for trt in qtl_keys: # for each trait
            trt_indices = self.all_indices(trt,list(data[:,0])) # indices of rows with this trait
            chrom = []
            pos = []
            val = []
            names = []
            for i in trt_indices: # for each row with this trait
                chrom.append(data[i,1]) # add chromosome information
                pos.append(data[i,2]) # add locus position information                
                val.append(data[i,3])
                names.append(data[i,4]) # add locus name information
            dic[trt] = [chrom,pos,val,names]
        return(dic)
        
    def pleiotropyManagement(self,data):
        """ This function manages pleitropy, i.e multi-traits QTL cases. We supose that there isn't the same position for two loci on the same chromosome for the same trait. """
        dic = self.dictGenerator(data)
        chrom = list(data[:,1])
        unique = np.unique(chrom,return_counts=True)
        unique_chrom = list(unique[0]) # unique values of chromosomes (sorted)
        chrom_count = list(unique[1]) # occurence of each chromosomes (sorted)
        
        for i in range(np.shape(data)[0]): # transformation of trait information from string to a list of string
            data[i,0] = [data[i,0]]
            data[i,3] = [data[i,3]]
            data[i,5] = [data[i,5]]
        
        for i in range(len(unique_chrom)): # for each unique value of chromosome
            if chrom_count[i] > 1: # if the chromosome is present more than once in the list
                polychrom = unique_chrom[i] # number of the chromosome
                print("polychrom",polychrom)
                chrom_index = self.all_indices(polychrom,list(data[:,1])) # indices of all the occurences of the chromosome
                print("chrom_index",chrom_index)
                polychrom_index = chrom_index[0] # index of the first occurence of the chromosome in the data
                size = len(chrom_index)-1 # number of loci positions to test for this chromosome
                j = 1 # because there is the first chromosome
                while j <= size :
                    print("j",j)
                    jj = chrom_index[j] # index of the row in the data for which the chromosome is also present
                    print("jj",jj)
                    if data[polychrom_index,2] == data[jj,2]: # if the position is the same (multi-trait QTL)
                        newname = data[polychrom_index,4]+","+data[jj,4] #
                        print("newname",newname)
                        for trt1 in data[polychrom_index,0]: # update of the dictionnary list
                            print("trt1",trt1)
                            for trt2 in data[jj,0]:
                                print("trt2",trt2)
                                pos1 = self.all_indices(polychrom,dic[trt1][0])[0]
                                print(pos1)
                                pos2 = self.all_indices(polychrom, dic[trt2][0])[0]
                                print(pos2)
                                print(dic[trt1][3])
                                print(dic[trt2][3])
                                dic[trt1][3][pos1] = newname
                                dic[trt2][3][pos2] = newname
                        data[polychrom_index,4] = data[polychrom_index,4]+","+data[jj,4]  # change the name of the first occurence in the data
                        data[polychrom_index,0].append(data[jj,0][0]) # add second trait in the list of trait of the first occurence
                        data[polychrom_index,3].append(data[jj,3][0])
                        data[polychrom_index,5].append(data[jj,5][0])
                        data = np.delete(data,jj,0) # delete the second occurence from data
                        chrom_index = list(np.add(chrom_index,-1)) # -1 as we deleted a row but may want to access another row
                    j = j+1
        return([data,dic])
    

