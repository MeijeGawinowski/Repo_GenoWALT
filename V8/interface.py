N_geno=int(input("How many genotypes do you want to generate ? "))
N_trait=int(input("How many traits do you want to study ? "))

list_trt=[]
for trt in range(N_trait):
    name_trt=input("Trait nb "+str(trt+1)+" : ")
    list_trt.append(name_trt) 

print(list_trt)

Reference_values={}
for i in range(N_geno):
    # for each genotype
    name_geno="Geno"+str(i)
    Reference_values[name_geno]={}
    for trt in list_trt :
        val=input("Value of "+trt+" for genotype nb "+str(i)+" : ")
        Reference_values[name_geno][trt]=val
        
print(Reference_values)