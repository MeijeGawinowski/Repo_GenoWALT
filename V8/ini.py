import pandas as pd
import numpy as np
import os
import random
import matplotlib.pyplot as plt
import csv

""" Warning, if you use this script you have to make sure that the traits of interest are the same in this file (for pheno_limits) than the parameters you will enter. If you are using WALTer you want to make sure that you have the right parameters in the variable and fixed genotypic parameters files. Also pay attention to the data QTL file you're using for the traits of interest """


pheno_limits = {'Param_PlHeight':[50,170],'L_B_max':[5,35],'GAI_c':[0,40]}
fileQTL = "fileQTL4.csv"


######################################################
### Functions that read QTL data from the QTL file ###
######################################################


def getDataQTL(fileQTL):
    file = open(fileQTL,"rb")
    df=pd.read_csv(file,delimiter=",",index_col=False)
    file.close()
    data = np.array(df)
    return(data)
    
def dictGenerator(data):
    """ This function generates a dictionnary list with each trait as key and list of chromosomes, list of loci positions and list of loci names as values. Pleiotropy isn't taken into account. """
    qtl_keys = list(np.unique(data[:,0])) # unique values of traits
    dic = {}
    for trt in qtl_keys: # for each trait
        trt_indices = all_indices(trt,list(data[:,0])) # indices of rows with this trait
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
    
def pleiotropyManagement(data):
    """ This function manages pleitropy, i.e multi-traits QTL cases. We supose that there isn't the same position for two loci on the same chromosome for the same trait. """
    dic = dictGenerator(data)
    chrom = list(data[:,1])
    unique = np.unique(chrom,return_counts=True)
    unique_chrom = list(unique[0]) # unique values of chromosomes (sorted)
    chrom_count = list(unique[1]) # occurence of each chromosomes (sorted)
        
    for i in range(np.shape(data)[0]): # transformation of trait information from string to a list of string
        if isinstance(data[i,0],list) == False :
            data[i,0] = [data[i,0]]
        if isinstance(data[i,3],list) == False :
            data[i,3] = [data[i,3]]
        if isinstance(data[i,5],list) == False :
            data[i,5] = [data[i,5]]
        
    for i in range(len(unique_chrom)): # for each unique value of chromosome
        if chrom_count[i] > 1: # if the chromosome is present more than once in the list
            polychrom = unique_chrom[i] # number of the chromosome
                #print("polychrom",polychrom)
            chrom_index = all_indices(polychrom,list(data[:,1])) # indices of all the occurences of the chromosome
                #print("chrom_index",chrom_index)
            polychrom_index = chrom_index[0] # index of the first occurence of the chromosome in the data
            size = len(chrom_index)-1 # number of loci positions to test for this chromosome
            j = 1 # because there is the first chromosome
            while j <= size :
                    #print("j",j)
                jj = chrom_index[j] # index of the row in the data for which the chromosome is also present
                    #print("jj",jj)
                if data[polychrom_index,2] == data[jj,2]: # if the position is the same (multi-trait QTL)
                    newname = data[polychrom_index,4]+","+data[jj,4] #
                        #print("newname",newname)
                    for trt1 in data[polychrom_index,0]: # update of the dictionnary list
                            #print("trt1",trt1)
                        for trt2 in data[jj,0]:
                                #print("trt2",trt2)
                            pos1 = all_indices(polychrom,dic[trt1][0])[0]
                                #print(pos1)
                            pos2 = all_indices(polychrom, dic[trt2][0])[0]
                                #print(pos2)
                                #print(dic[trt1][3])
                                #print(dic[trt2][3])
                            dic[trt1][3][pos1] = newname
                            dic[trt2][3][pos2] = newname
                    data[polychrom_index,4] = data[polychrom_index,4]+","+data[jj,4]  # change the name of the first occurence in the data
                    data[polychrom_index,0].append(data[jj,0][0]) # add second trait in the list of trait of the first occurence
                    data[polychrom_index,3].append(data[jj,3][0])
                    data[polychrom_index,5].append(data[jj,5][0])
                    data = np.delete(data,jj,0) # delete the second occurence from data
                    chrom_index = list(np.add(chrom_index,-1)) # -1 as we deleted a row but may want to access another row
                j = j+1
    return(data)
    
def all_indices(value,qlist):
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
        
def limits_definer(dict):
    list_trt=list(dict.keys())
    new_dict={}
    for i in range(len(list_trt)):
        trt=list_trt[i]
        #print("trt=",trt)
        values=dict[trt][2]
        #print("values=",values)
        neg=[]
        pos=[]
        for j in range(len(values)):
            if values[j] < 0:
                neg.append(values[j])
            if values[j] > 0:
                pos.append(values[j])
            #print(neg,pos)
        min_trt=sum(neg)
        max_trt=sum(pos)
        new_dict[trt]=[min_trt,max_trt]
    return(new_dict)    

def conversion_allele_QTL(data): 
# data is the QTL table
    N_qtl = np.shape(data)[0]
    dico = {}
    for qtl in range(N_qtl): # for each qtl in data_qtl
        values_qtl = data[qtl,3]
        # for each type qtl
        subdict = {}
        subdict["[0,0]"] = []
        subdict["[0,1]"] = []
        subdict["[1,0]"] = []
        subdict["[1,1]"] = []
        #print(values_qtl)
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
        dico[qtl] = subdict
    return(dico)  




#############################################
### Functions that for genotypes sampling ###
#############################################


def random_product(*args, **kwds):
    "Random selection from itertools.product(*args, **kwds)"
    pools = args * kwds.get('repeat', 1)
    return([random.choice(pool) for pool in pools])
    
    
def samplingGeno(samp_size,N_qtl,heterozygoty):
    
    if heterozygoty == 0: # only homozygotes
        alleles_list = ["[0,0]","[1,1]"]
    if heterozygoty == 1: # homozygotes and heterozygotes
        alleles_list = ["[0,0]","[0,1]","[1,0]","[1,1]"]
    if heterozygoty == 2: # only heterozygotes
        alleles_list = ["[0,1]","[1,0]"]
    
    stock=[]
    i=0
    while i <= samp_size :    
        res=random_product(alleles_list,repeat=N_qtl)
        dico={}
        for e in range(len(res)):
            dico[e]=res[e]
        stock.append(dico)
        i=i+1
    return(stock)
    


###################################################
### Functions for genotype/phenotype conversion ###
###################################################

def getListQTL(dico_geno,conv):
    N_qtl = len(dico_geno)
    list_out = []
    for j in range(N_qtl):
        idx=list(dico_geno.keys())[j]
        allel = dico_geno[idx] # allel value for the individual i at the qtl j       
        al1=str(allel).replace(" ","")
        al2="["+str(allel[1])+","+str(allel[4])+"]"
        # conv[j] est liste allèles/valeurs pour le qtl j
        list_allel = list(conv[j].keys()) # alleles list for qtl j
        list_valqtl = list(conv[j].values()) # qtl j values list for each allele
        try :
            idx = list_allel.index(al1) # index of the allele in the list
        except ValueError :
            idx = list_allel.index(str(list(allel)))
        val = list_valqtl[idx] # qtl value matching this allele
        list_out.append(val)
    return(list_out)

def offspringQTL(dictGeno,conv):
    N_qtl = len(dictGeno)
    list_out = []
    for j in range(N_qtl):
        idx=list(dictGeno.keys())[j]
        allel = dictGeno[idx] # allel value for the individual i at the qtl
        al1=str(allel).replace(" ","")       
        al2="["+str(allel[1])+","+str(allel[4])+"]"
            # conv[j] est liste allèles/valeurs pour le qtl j
        list_allel = list(conv[j].keys()) # alleles list for qtl j
        list_valqtl = list(conv[j].values()) # qtl j values list for each allele
        try :
            idx = list_allel.index(al1) # index of the allele in the list
        except ValueError :
            idx = list_allel.index(str(list(allel)))
        val = list_valqtl[idx] # qtl value matching this allele
        list_out.append(val)
    return(list_out)
        
def rescaling(x,old_lims,new_lims):
    old_min = old_lims[0]
    old_max = old_lims[1]
    new_min = new_lims[0]    
    new_max = new_lims[1]
    y=((new_min-new_max)/(old_min-old_max))*(x-old_min)+new_min
    return(y)
    
def trt_indices(tab_qtl,trt):
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


#######################
### QTL data import ###
#######################

# import QTL data
dataQTL1 = getDataQTL(fileQTL)
#print(dataQTL1)
dict_QTL=dictGenerator(dataQTL1)
dataQTL2 = pleiotropyManagement(dataQTL1)
#print(dataQTL2)
N_qtl=np.shape(dataQTL2)[0]
# set conversion between alleles and values for each QTL
conversion_dico = conversion_allele_QTL(dataQTL2)
qtl_limits=limits_definer(dict_QTL)


################################################################
### Choice of genotypes by the user (N_trait, N_ind, values) ###
################################################################

N_geno=int(input("How many genotypes do you want to generate ? "))
N_trait=int(input("How many traits do you want to study ? "))

list_trt=[]
for trt in range(N_trait):
    name_trt=input("Trait nb "+str(trt+1)+" : ")
    list_trt.append(name_trt) 

#print(list_trt)

Reference_values={}
names_geno=[]
N_ind=0
new_dict={}
for i in range(N_geno):
    # for each genotype
    name_geno="Geno"+str(i)
    new_dict[name_geno]={}
    names_geno.append(name_geno)
    Reference_values[name_geno]={}
    for trt in list_trt :
        val=float(input("Value of "+trt+" for genotype nb "+str(i)+" : "))
        Reference_values[name_geno][trt]=val
    prop=input("How many individuals for this genotype ? ")
    heterozygoty_value=int(input("What type of loci do you want ? \n 0: Only homozygotes \n 1: Homozygotes and heterozygotes \n 2: Only heterozygotes \n"))
    N_ind+=int(prop)
    new_dict[name_geno]['eff']=int(prop)
    new_dict[name_geno]['heterozygoty'] = heterozygoty_value
    
#print("new_dict=",new_dict)

N_samp=int(input("How many samples for genotype search ? "))



###############################################################
### Search of best genotypes and writing of geno_conversion ###
###############################################################

# initialization of geno_conversion
geno_conversion=np.empty([1,2],dtype=object)

# computation of best genotypes for each wanted phenotype
#print("NAMES=",names_geno)
for geno_name in names_geno:
    #print("geno=",geno_name)
    input_values = Reference_values[geno_name]
    chosen_heterozygoty = new_dict[geno_name]['heterozygoty']
    # sample genotypes among all possibilities
    sampled_geno=samplingGeno(N_samp,N_qtl,chosen_heterozygoty)
    #print(sampled_geno)
    
    list_interest=list(input_values.keys())
    wanted_values=list(input_values.values())
    
    
    # initialization of the dictionnary list of results
    sampled_results={}
    for i in range(len(sampled_geno)):
        sampled_results[i]={}
        for trt in list_interest:
            sampled_results[i][trt]=0
    
    #print("First sampled results are ",sampled_results)
    for trt in list_interest:
            #print("current trait is ",trt)
            trt_idx=list_interest.index(trt)
            list_trt_indices=trt_indices(dataQTL2,trt)
            #print("indices of values for this trait are : ",list_trt_indices)
            for i in range(len(sampled_geno)):
                geno=sampled_geno[i]
                list_qtl=getListQTL(geno,conversion_dico)
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
                sum_qtl=sum(list_qtl)
                #print(sum_qtl,qtl_limits[trt],pheno_limits[trt])
                val_trt=rescaling(sum_qtl,qtl_limits[trt],pheno_limits[trt])
                sampled_results[i][trt]=val_trt
                    
    #print("sampled results are",sampled_results)
    
    diff_list=[]
    for i in range(len(sampled_geno)):
        subdict=sampled_results[i]
        tested_values=list(subdict.values())
        diff_val=[abs(ref_elt-test_elt) for ref_elt, test_elt in zip(wanted_values,tested_values) ]
        diff_val=sum(diff_val)
        diff_list.append(diff_val)
        
    #print("diff_list=",diff_list)
    
    min_elt=min(diff_list)
    print("For genotype "+geno_name+" the minimum difference is of ",min_elt)
    min_idx=diff_list.index(min_elt)
    chosen_geno=sampled_geno[min_idx]
    #print("chosen genotype is ",chosen_geno)
    chosen_val=sampled_results[min_idx]
    new_dict[geno_name].update(chosen_val)
    new_dict[geno_name]['dico']=chosen_geno
    #print("Now new_dict=",new_dict)
    print("For genotype "+geno_name+" the closest value is of ",chosen_val)
    
    new_raw=np.array([geno_name,chosen_geno])
    #print(new_raw)
    geno_conversion=np.row_stack([geno_conversion,new_raw.astype(np.object)])

geno_conversion=np.delete(geno_conversion,0,0)
#print(geno_conversion)

geno_file = open("first_conversion.csv","w")
geno_writer = csv.writer(geno_file,delimiter="\t")
geno_writer.writerows(geno_conversion)
geno_file.close()

##################################
### Writing of population_data ###
##################################

# writing of population data
final_data=np.empty((N_ind+1,N_trait+3),dtype=object)
#print("SHAPE=",np.shape(final_data))
# filling final_data header
final_data[0,0]='individual'
final_data[0,1]='fitness'
final_data[0,2]='genotype'
for trt in range(N_trait):
    #print("TRT=",trt)
    final_data[0,trt+3]=list_interest[trt]
        
    i=1 # row 1 because 0 is header
    for k in range(len(names_geno)):
        geno_name=names_geno[k]
        #print("GENO=",geno_name)
        n_ind = new_dict[geno_name]['eff']
        #print("n_ind=",n_ind)
        lim=n_ind+i-1
        while i <= lim:
            #print("i=",i)
            final_data[i,0]=i-1 # index of the individual in the population starting at 0
            final_data[i,2]=geno_name
            for trt in range(N_trait):
                trt_name=list_interest[trt]
                final_data[i,trt+3]=new_dict[geno_name][trt_name]
                #print("val=",new_dict[geno_name][trt_name])    
                # calculation of fitness
                if trt_name=="Param_PlHeight":
                    final_data[i,1]=new_dict[geno_name][trt_name]
            i=i+1

        
        #print("ind=",ind,len(ind))
        #print("final geno=",final_geno,len(final_geno))
        #print("final fitness=", final_fitness,len(final_fitness))
        
        #final_data = np.column_stack([ind,list_final_geno,np.array(final_fitness).astype(np.object)]) 
        
#print("final data is : ",final_data)
          
# Writing new table of genotype names
pop_file = open("first_population.csv","w")
pop_writer = csv.writer(pop_file,delimiter="\t")
pop_writer.writerows(final_data)
pop_file.close()
