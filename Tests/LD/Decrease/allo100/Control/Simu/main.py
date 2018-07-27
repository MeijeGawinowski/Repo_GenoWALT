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
from simuPOP.utils import export, saveCSV

os.chdir('/home/meije/Documents/LD_Tests/Decrease/Control/Simu')

from Reader import *
from WALTerReader import *
from Population import *
from Simulation import *
from Offspring import *


# Import of QTL data
qtl_object = Reader()
qtl_data = qtl_object.dataQTL()
tab_qtl = qtl_data[0][0] # array with the genetic positions of QTLs
dict_qtl = qtl_data[0][1] # dictionary list of the QTLs
tab_loci = qtl_data[1] # genetic positions of other loci
# print("QTL table : ", tab_qtl)
# print("QTL dictionnary : ", dict_qtl)
# print("Loci table : ", tab_loci)



# Import of WALTer data (genotype and phenotype)
init_walter_object = WALTerReader("Dict_LD1.txt")
init_data_walter = init_walter_object.dataPheno(tab_qtl)
N_ind = init_data_walter[0] # number of individuals in the population
init_genoPop = init_data_walter[1] # population genotype
conv = init_data_walter[2] # conversion table between alleles and QTL values




# Creation of the parent population 
init_pop_object = Population()
init_pop = init_pop_object.Creation(tab_qtl,tab_loci,init_genoPop,N_ind)
print("Generation nb 0")
sim.dump(init_pop)

N_gen = 25

parpop = init_pop # initialization of the parental population (gen 0)
saveCSV(init_pop,sexFormatter=None,affectionFormatter=None,genoFormatter={(0,0):"1|1",(1,0):"2|1",(0,1):"1|2",(1,1):"2|2"},sep="\t",filename="/home/meije/Documents/LD_Tests/Decrease/Control/thetaD/Rep5/pop0.txt")
for g in range(1,N_gen):
	print("Generation nb ",g)
	# Reproduction of the parent population
	simu_object = Simulation(N_ind)
	offpop = simu_object.create_controlledSimulation(parpop)
	namefile = "/home/meije/Documents/LD_Tests/Decrease/Control/thetaD/Rep5/pop"+str(g)+".txt"
	saveCSV(offpop,sexFormatter=None,affectionFormatter=None,genoFormatter={(0,0):"1|1",(1,0):"2|1",(0,1):"1|2",(1,1):"2|2"},sep="\t",filename=namefile)
	off_object = Offspring("")
	geno_off = off_object.Result(offpop,conv,tab_qtl) # dictionary of new pop
	parpop_object = Population()
	parpop = parpop_object.Creation(tab_qtl,tab_loci,geno_off,N_ind)
	# sim.dump(parpop)
