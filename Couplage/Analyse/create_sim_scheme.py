import os
import csv
import numpy as np

os.chdir('/home/meije/Documents/GenoWALT/Couplage')

# ouverture fichier avec parametres genotypiques fixes
file = open('fixed_geno_params.csv',"rb")
df=pd.read_csv(file,delimiter="\t")
file.close()
fixed_geno_params_dict=df.to_dict(orient='records')[0]

# ouverture fichier avec parametres genotypiques variables
file2=open("geno_ini.csv","rb")
df2=pd.read_csv(file2,delimiter="\t",index_col=False)
file2.close()
variable_geno_params_dict=df2.to_dict(orient='records')[0]

# ouverture fichier avec parametres experimentaux fixes
file3=open("fixed_expe_cond_params.csv","rb")
df3=pd.read_csv(file3,delimiter="\t",index_col=False)
file3.close()
fixed_expe_cond_params_dict=df3.to_dict(orient='records')[0]



# creation liste dictionnaire pour parametres fixes genotype
string_geno=variable_geno_params_dict["genotypes"]
list_geno=string_geno.split(",")
fixed_params = list(fixed_geno_params_dict.keys())
geno_params_dict={}
for geno in list_geno:
    print(geno)
    for par in fixed_params:
        name=par+"_"+geno
        val=fixed_geno_params_dict[par]
        geno_params_dict[name]=val

# mise a jour du dico avec les trois fichiers
geno_params_dict.update(variable_geno_params_dict)
geno_params_dict.update(fixed_expe_cond_params_dict)

# recuperation keys et values du dico
keys=list(geno_params_dict.keys())
values=list(geno_params_dict.values())

# ecriture dans un fichier
scheme_tab=np.array([keys,values])
scheme_file = open("sim_scheme.csv","w")
scheme_writer = csv.writer(scheme_file,delimiter="\t")
scheme_writer.writerows(scheme_tab)
scheme_file.close()



