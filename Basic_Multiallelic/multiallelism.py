import itertools as iter
import math as math

distrib = "expo"

dat=dataqtl[0][0]
N_qtl=np.shape(dat)[0]
values_qtl=dat[:,3]
nb_allele = dat[:,6]

def affine(x,N,v):
    return((x*v)/(N-1))
    
def expo(x,N,v):
    return((v*(math.exp(x)-1))/(math.exp(N-1)-1))


dict={}
for qtl in range(N_qtl):
    print("qtl ",qtl)
    subdict = {}
    nb=nb_allele[qtl]
    print(nb)
    v=values_qtl[qtl]
    print("v = ",v)
    if distrib == "affine":
        values=[affine(i,nb,v) for i in range(nb)] # allele values for qtl
        print("values = ",values)
    if distrib == "expo":
        values=[expo(i,nb,v) for i in range(nb)]
        print("values=",values)
    arr=list(iter.product(list(range(nb)),repeat=2))
    for k in range(len(arr)):
        print("k=",k)
        arr_elt = list(arr[k]) # list(tuple) avec les deux allèles
        al1 = arr_elt[0]
        al2 = arr_elt[1]
        meanval = (values[al1]+values[al2])/2
        print("res : ",al1,al2,meanval)
        subdict[str(arr_elt)] = meanval
    dict[qtl]=subdict
   
   

