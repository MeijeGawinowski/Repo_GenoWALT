import numpy as np
import csv

def pheno2sim_scheme(data_pheno,list_trt):
    N_ind=np.shape(data_pheno)[0]
    Size=np.shape(data_pheno)[1]
    print("size=",Size)
    geno_col=data_pheno[:,2][1:N_ind]
    print(geno_col)
    N_trait=Size-3
    print(" and Nb trait is ",N_trait)
    tab_unique=np.unique(geno_col,return_index=True,return_counts=True)
    print(tab_unique)
    names=list(tab_unique[0])
    N_geno=len(names)
    string_names=",".join(names)
    print(string_names)
    indices=list(tab_unique[1])
    prop=list(tab_unique[2])
    print(prop)
    string_prop=",".join(str(x) for x in prop)
    print(string_prop)
    list_val_trt=[]

    final_scheme=np.empty([2,2],dtype=object)
    final_scheme[0,0]="genotypes"
    final_scheme[0,1]="genotype_numbers"
    final_scheme[1,0]= string_names
    final_scheme[1,1]= string_prop
    

    
    for gidx in range(len(names)):
        print("gidx=",gidx)
        # for each genotype
        for trt in range(N_trait):
            name_trt=data_pheno[0,trt+3]
            colname=name_trt+"_"+names[gidx]
            idx=indices[gidx]
            val_trt=data_pheno[idx+1,trt+3]

            new_col=np.array([colname,val_trt],dtype=object)
            final_scheme=np.column_stack([final_scheme,new_col.astype(np.object)])
            print("final_scheme")
            #final_scheme[0,2+gidx+trt]=colname
            #final_scheme[1,2+gidx+trt]=val
    
    
    
    # Writing of a scheme file in a csv
    scheme_file = open("first_variable_geno_params.csv","w")
    scheme_writer = csv.writer(scheme_file,delimiter="\t")
    scheme_writer.writerows(final_scheme)
    scheme_file.close()
        
    return(final_scheme)

    
        
