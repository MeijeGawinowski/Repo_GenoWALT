#
# main.py $
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

import os

#os.chdir('/home/meije/Documents/GenoWALT/V7_PH')

from Reader import *
from WALTerReader import *
from Population import *
from Simulation import *
from Offspring import *
from pheno_sim import *
from simuPOP.utils import export, saveCSV


# Import of QTL data
qtl_object = Reader()
qtl_data = qtl_object.dataQTL()
tab_qtl = qtl_data[0] # array with the genetic positions of QTLs
dict_qtl = qtl_data[1] # dictionary list of the QTLs
limits_qtl_dict=qtl_object.limits_definer(dict_qtl)
print("QTL table : ", tab_qtl)
print("QTL dictionnary : ", dict_qtl)
print("Dictionnary of QTL values limits : ",limits_qtl_dict)






 

N_gen=1
par_WALTer=True

for g in range(N_gen):
    print("g=",g)
    # Import of WALTer data (genotype and phenotype)
    walter_object = WALTerReader()
    data_walter = walter_object.dataWALTer(tab_qtl)
    genoPop = data_walter[0] # population genotype
    print("Parent population genotype : ", genoPop)
    conv = data_walter[1] # conversion table between alleles and QTL values
    print("Conversion table QTL/allel : ", conv)
    qtl_tab = data_walter[2]
    print("Table des QTL : ",qtl_tab)
    header_list=data_walter[3]
    print("Header : ",header_list)
    fitness_list=data_walter[4]
    print("Fitness : ",fitness_list)
    walter_array=data_walter[5]
    print("Data WALTer : ",walter_array)

    # Creation of the parent population
    parpop_object = Population()
    parpop = parpop_object.Creation(tab_qtl,genoPop,fitness_list,header_list,walter_array)
    #print("pop n",0," is :")
    #sim.dump(parpop)
    print("Info fields are ",list(parpop.infoFields()),list(parpop.indInfo('fitness')))
    print("pop is")
    sim.dump(parpop)
    # Reproduction of the parent population
    simu_object = Simulation()
    offpop = simu_object.create_FitnessSimulation(parpop)
    off_object = Offspring(qtl_object.list_trt,limits_qtl_dict)
    
    pheno_data = off_object.GenoWALT2WALTer(offpop,tab_qtl,conv)
    print("pheno_data is ",pheno_data)
    
    if par_WALTer == True:
        pheno2sim_scheme(pheno_data,qtl_object.list_trt)




