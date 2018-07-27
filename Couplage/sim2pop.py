import csv
import pandas as pd
import numpy as np


# open sim_scheme file

file = open("first_variable_geno_params.csv","rb")
df=pd.read_csv(file,delimiter="\t",index_col=False)
file.close()
scheme_header=list(df)
scheme_values = np.array(df)
scheme_values = list(scheme_values[0,:])
print("scheme header is ",scheme_header)
print("scheme values are ",scheme_values)

idx_genotypes = scheme_header.index('genotypes')
str_genotypes = scheme_values[idx_genotypes]
list_genotypes = str_genotypes.split(",")
print("List of genotypes is ",list_genotypes)

# open population_data file
file = open("first_population.csv","rb")
df=pd.read_csv(file,delimiter="\t",index_col=False)
file.close()
pop_header=list(df)
pop_array=np.array(df)
list_geno_pop = pop_array[:,2]
print("List of genotypes in the population is ",list_geno_pop)


list_interest = ['Param_PlHeight','L_B_max','GAI_c']


for trt in list_interest :
	print("TRT=",trt)
	new_col = []
	for name_geno in list_geno_pop:
		print("GENO=",name_geno)
		name_param = trt+"_"+name_geno
		param_idx = scheme_header.idx(name_param)
		param_value = scheme_values[param_idx]
	print("new col is ",new_col, len(new_col))
	pop_header.append(trt)
	pop_array = np.column_vstack([pop_array,np.array(new_col).astype(np.object)])
	
print("Now pop_array is ",pop_array)
print("header=",pop_header)

