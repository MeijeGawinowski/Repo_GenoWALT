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


    def __init__(self,trt_interest,qtl_limits):
        self.trt_interest = trt_interest
        self.qtl_limits=qtl_limits
        self.pheno_limits = {'Param_PlHeight':[50,170],'L_B_max':[5,35],'GAI_c':[0,40]}
        self.trt_mean = {'Param_PlHeight': 76.1}
        self.choice_pheno = "rescaling" # "addition" or "rescaling"


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
        #print("Genotab",tab_out)
        return(tab_out)
        
    def dictGenoOffspring(self,genotab):
        dict = {}
        for i in range(np.shape(genotab)[0]):
            subdict = {}
            for j in range(np.shape(genotab)[1]):
                subdict[j] = str(genotab[i,j])
            dict[i] = subdict
        #print("Offspring geno dict ",dict)
        return(dict)

    def offspringQTL(self,dictGeno,conv):
        N_qtl = len(dictGeno)
        list_out = []
        for j in range(N_qtl):
            idx=list(dictGeno.keys())[j]
            allel = dictGeno[idx] # allel value for the individual i at the qtl
            al1=str(allel).replace(" ","")       
            al2="["+str(allel[1])+","+str(allel[4])+"]"
            # conv[j] est liste alleles/valeurs pour le qtl j
            list_allel = list(conv[j].keys()) # alleles list for qtl j
            list_valqtl = list(conv[j].values()) # qtl j values list for each allele
            try :
                idx = list_allel.index(al1) # index of the allele in the list
            except ValueError :
                idx = list_allel.index(str(list(allel)))
            val = list_valqtl[idx] # qtl value matching this allele
            list_out.append(val)
        return(list_out)
        
    def trt_indices(self,tab_qtl,trt):
        sublist = []
        list_idx = []
        for i in range(np.shape(tab_qtl)[0]):
            list_trt = tab_qtl[i,0]
            try :
                idx = list_trt.index(trt)
                sublist.append(idx)
                list_idx.append(i)
            except ValueError :
                sublist.append(-1)
        return(sublist)
    
    def rescaling(self,x,old_lims,new_lims):
        old_min = old_lims[0]
        old_max = old_lims[1]
        new_min = new_lims[0]
        new_max = new_lims[1]
        y=((new_min-new_max)/(old_min-old_max))*(x-old_min)+new_min
        return(y)
        
    
    def offspringPheno(self,list_qtl,trt):
        if self.choice_pheno == "addition" :
            val_trt = self.trt_mean[trt] + sum(list_qtl)
        if self.choice_pheno == "rescaling":
            val_trt = sum(list_qtl)
            val_trt = self.rescaling(val_trt,self.qtl_limits[trt],self.pheno_limits[trt])
        return(val_trt)
        
        
                
    def GenoWALT2WALTer(self,offpop,tab_qtl,conv):
        genotab=self.offspringGeno(offpop,conv,tab_qtl)
        res=self.dictGenoOffspring(genotab)
        N_ind=len(res)
        N_trait = len(self.trt_interest)
                
        test=[False]*N_ind
        #print("test=",test)
        new_dict={}
        for i in range(N_ind):
            if test[i]==False:
                same_geno=[i]
                test[i]=True
                for j in range(N_ind):
                    if i != j:
                        if res[i] == res[j]:
                            same_geno.append(j)
                            test[j]=True
                name="Geno"+str(i) 
                eff=len(same_geno)
                dict_geno=res[i]
                new_dict[name]={}
                new_dict[name]['eff']=eff
                new_dict[name]['dico']=dict_geno
        
        print("new_dict=",new_dict) 
        N_geno=len(new_dict)
        names=list(new_dict.keys())
        #print(N_geno,names)
        string_names=",".join(names)
        #print("str names=",string_names)
        tot_names=[]
        
        string_prop=[]
        int_prop=[]
        for n in range(len(names)):
            key=names[n]
            n_name=[key]*N_trait
            tot_names=tot_names+n_name
            eff=new_dict[key]['eff']
            int_prop.append(eff)
            string_prop.append(str(eff))
        #print("prop=",prop)
        string_prop=",".join(string_prop)
        #print("str prop=",string_prop)
        
        final_geno=np.empty([N_geno,2],dtype=object)
                    
        # Remplir final_geno
        for g in range(N_geno):
            key=names[g]
            final_geno[g,0]=key
            final_geno[g,1]=new_dict[key]['dico']

        # Calculer les phenotypes
        
        for trt in self.trt_interest:
            print("current trait is ",trt)
            trt_idx=self.trt_interest.index(trt)
            list_trt_indices=self.trt_indices(tab_qtl,trt)
            print("indices of values for this trait are : ",list_trt_indices)
            for geno_name in names:
                dict_geno=new_dict[geno_name]['dico']
                list_qtl=self.offspringQTL(dict_geno,conv)
                N_qtl=len(list_qtl)
                comp=0
                qtl=0
                while comp <= N_qtl-1:
                    if list_trt_indices[comp] == -1: # no QTL for trt
                        del list_qtl[qtl]
                        comp = comp+1
                    else:
                        qtl_idx = list_trt_indices[qtl]
                        list_qtl[qtl]=list_qtl[qtl][qtl_idx]
                        qtl = qtl+1
                        comp=comp+1
                val_trt=self.offspringPheno(list_qtl,trt)
                new_dict[geno_name][trt]=val_trt
                
        print("Now new dict is : ",new_dict)
        
        
        
        N_ind=sum(int_prop)
        print("Nb of individuals is ",N_ind)
        final_data=np.empty((N_ind+1,N_trait+3),dtype=object)
        
        # filling final_data header
        final_data[0,0]='individual'
        final_data[0,1]='fitness'
        final_data[0,2]='genotype'
        for trt in range(N_trait):
            final_data[0,trt+3]=self.trt_interest[trt]
        
        i=1 # row 1 because 0 is header
        for geno_name in names:
            n_ind = new_dict[geno_name]['eff']
            lim=n_ind+i-1
            print("eff is ",n_ind," and lim is ",lim)
            while i <= lim:
                print("i=",i)
                final_data[i,0]=i-1 # index of the individual in the population starting at 0
                final_data[i,2]=geno_name
                for trt in range(N_trait):
                    trt_name=self.trt_interest[trt]
                    final_data[i,trt+3]=new_dict[geno_name][trt_name]
                    
                    # calculation of fitness
                    if trt_name=="Param_PlHeight":
                        final_data[i,1]=new_dict[geno_name][trt_name]
                i=i+1

        print("final data is : ",final_data)
       
        # Writing of a scheme file in a csv
        #data_file = open("population_data.csv","w")
        #data_writer = csv.writer(data_file,delimiter="\t")
        #data_writer.writerow(['individual','fitness','genotype'])
        #data_writer.writerows(final_data)
        #data_file.close()
        
        # Writing new table of genotype names
        geno_file = open("geno_conversion.csv","w")
        geno_writer = csv.writer(geno_file,delimiter="\t")
        geno_writer.writerows(final_geno)
        geno_file.close()
        
        return(final_data)
    
        

    
