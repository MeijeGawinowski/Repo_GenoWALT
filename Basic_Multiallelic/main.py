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

os.chdir('/home/meije/Documents/V5')

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
print("QTL table : ", tab_qtl)
print("QTL dictionnary : ", dict_qtl)
print("Loci table : ", tab_loci)



# Import of WALTer data (genotype and phenotype)
walter_object = WALTerReader()
data_walter = walter_object.dataPheno(tab_qtl)
genoPop = data_walter[0] # population genotype
print("Parent population genotype : ", genoPop)
conv = data_walter[1] # conversion table between alleles and QTL values
print("Conversion table QTL/allel : ", conv)
pheno = data_walter[2] # phenotype object with the population phenotype table, the header (list of the different traits) and the number of individuals in the population
print("Parent population phenotype : ", pheno)



# Creation of the parent population 
parpop_object = Population()
parpop = parpop_object.Creation(tab_qtl,tab_loci,dict_qtl,genoPop,pheno)
mean_PH = parpop_object.getMean(parpop)
print("Parent population :")
sim.dump(parpop)
print("PH mean = ", mean_PH)
print("Information Fields of the parent population : ",parpop.infoFields())



# Reproduction of the parent population
simu_object = Simulation()
offpop = simu_object.create_controlledSimulation(parpop)
print("Offsprings population : ")
sim.dump(offpop)
print("Information Fields of the offspring population : ",offpop.infoFields())



# Retrieval of genotypic and phenotypic data of the offspring population
off_object = Offspring()
res = off_object.Result(offpop,conv,dict_qtl,mean_PH,qtl_object.choice)
print("Genotype dictionary of the offspring population (dictRes.txt file) : ", res[0])
print("Phenotype of the traits of interest of the offspring population (phenoRes.csv file): ", res[1])

modeTagger=offpop.indInfo('mode')
print("Reproduction mode : ",modeTagger)
fathers=offpop.indInfo('father_idx')
print("Fathers ids : ",fathers)
mothers=offpop.indInfo('mother_idx')
print("Mothers ids : ",mothers)
print("offpsring pop size = ",offpop.popSize())
