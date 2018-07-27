#
# Simulation.py $
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

import math as math
import numpy as np
import simuPOP as sim
import copy as cp


class Simulation():
    
    def __init__(self):
        self.generations = 10
        self.replicat = 1
        self.seed = 1
        self.selfing = 0.8 # proportion of autogamy in the population
        self.selection = 'fitness' # selection trait for reproduction
        self.opt = "Haldane" # computation option for recombination rates
        self.newSize = 10

            
    def create_FitnessSimulation(self,parpop):
        """ This function runs the simulation """
        self.initOps = self.create_InitOps()
        self.preOps = self.create_PreOps()
        self.mating = self.create_Mating(parpop)
        self.postOps = self.create_PostOps()
        self.finalOps = self.create_FinalOps()
        
		#initialization of new information fields
        parpop.addInfoFields(['mode','mother_idx','father_idx'])
        parpop.evolve(initOps=self.initOps, preOps=self.preOps, matingScheme=self.mating, finalOps=self.finalOps,gen=self.generations)

        return(parpop)
    


    def create_controlledSimulation(self,parpop):
        N_loc = sum(list(parpop.numLoci()))

        recomb_loci = range(N_loc-1) # for the computation of recombination rate per locus we don't take into account the last locus (because it has already been taken into account by the computation of previous loci)
        
        Des = list(parpop.indInfo('offspring'))
                
        N_ind = parpop.popSize()
        ind_idx=list(range(N_ind))
        
        DesTot=sum(Des)
        prop=[]
        statN_des=[]
        for i in range(len(Des)):
            p = Des[i]/DesTot
            prop.append(p)
            n = round(self.newSize*p)
            statN_des.append(n)
        print("Nb of offspring per parent : ",statN_des)
        dynN_des=cp.deepcopy(statN_des) # updated during reproduction process    
        newpop=cp.deepcopy(parpop)
        newpop.removeIndividuals(indexes=range(N_ind))
        newpop.removeInfoFields(list(parpop.infoFields()))
        newpop.addInfoFields(['mode','mother_idx','father_idx'])
        newpop.setSubPopName('test',0)
            
        for i in range(N_ind):
            # for the individual i
            print("Progeny of individual N°",i)
            scheme=list(np.random.choice([0,1],size=dynN_des[i], replace=True,p=[0.8,0.2]))
            print(scheme)
            N_auto=scheme.count(0)
            N_allo=scheme.count(1)
            print("allof=",N_allo," and autof=",N_auto)
                
            # production of offsprings by allofecondation
            print("ALLOF")
            for j in range(N_allo): # for the j offsprings produced by allofecondation by individual i
                offsp = cp.deepcopy(dynN_des)
                offsp[i] = 0 # sure not choose i as partner at the next step
                mates_idx = self.check_list(offsp) # choice of a partner with a non null offspring value
                if len(mates_idx) == 0: # there is no partner left
                    buf_pop=parpop.extractIndividuals(i)
                    buf_pop.removeInfoFields(list(parpop.infoFields()))
                    buf_pop.evolve(matingScheme=sim.SelfMating(ops = [sim.Recombinator(rates=self.calcRecombination(parpop,recomb_loci),loci=recomb_loci)],subPopSize=1),gen=1)
                    dynN_des[i]=dynN_des[i]-1
                    buf_pop.addInfoFields(['mode','mother_idx','father_idx'])
                    buf_pop.setIndInfo(0,'mode')
                    buf_pop.setIndInfo(i,'mother_idx')
                    buf_pop.setIndInfo(-1,'father_idx')
                    buf_pop.setSubPopName('test',0)
                    newpop.addIndFrom(buf_pop)
                else :
                    # recalculation of proportions
                    mates_prop=[prop[j] for j in mates_idx]
                    s=sum(mates_prop)
                    mates_prop=[mates_prop[i]/s for i in range(len(mates_prop))]
                        
                    # sampling of a partner within the ones possible ponderated by their proportions
                    lover=int(np.random.choice(mates_idx,size=1,replace=True,p=mates_prop))
                    print("Partner nb",j," for ",i," is ind nb",lover)
                    
                    dynN_des[lover]=dynN_des[lover]-1
                    dynN_des[i]=dynN_des[i]-1
                    # buffer population with 2 parent individuals
                    buf_pop=parpop.extractIndividuals([i,lover])
                    buf_pop.removeInfoFields(list(parpop.infoFields()))
                    buf_pop.evolve(matingScheme=sim.HermaphroditicMating(allowSelfing = False,ops=[sim.Recombinator(rates=self.calcRecombination(parpop,recomb_loci),loci=recomb_loci)],subPopSize=1),gen=1)
                    buf_pop.addInfoFields(['mode','mother_idx','father_idx'])
                    buf_pop.setIndInfo(1,'mode')
                    buf_pop.setIndInfo(i,'mother_idx')
                    buf_pop.setIndInfo(lover,'father_idx')
                    buf_pop.setSubPopName('test',0)
                    newpop.addIndFrom(buf_pop)
            
            
            print("AUTOF")
            # production of offsprings by autofecondation
            buf_pop=parpop.extractIndividuals(i) # parent population with one individual
            buf_pop.removeInfoFields(list(parpop.infoFields()))
            buf_pop.evolve(matingScheme=sim.SelfMating(ops = [sim.Recombinator(rates=self.calcRecombination(parpop,recomb_loci),loci=recomb_loci)],subPopSize=N_auto),gen=1)
            dynN_des[i]=dynN_des[i]-N_auto
            #print("N_des after trt of i",dynN_des)
            buf_pop.addInfoFields(['mode','mother_idx','father_idx'])
            buf_pop.setIndInfo([0]*N_auto,'mode')
            buf_pop.setIndInfo([i]*N_auto,'mother_idx')
            buf_pop.setIndInfo([-1]*N_auto,'father_idx')
            buf_pop.setSubPopName('test',0)
            newpop.addIndFrom(buf_pop)
            
        newpop.mergeSubPops()
        return(newpop)
        
        
    def create_InitOps(self):
        initOps = sim.InitSex()
        return(initOps)
    
    def create_PreOps(self):
        preOps = []
        return(preOps)
    
    def create_PostOps(self):
        postOps = []
        return(postOps)
    
    def create_FinalOps(self):
        finalOps = []
        return(finalOps)
        
    def alloTagger(self,mode):
        return(1)
    
    def autoTagger(self,mode):
        return(0)
    
    def check_list(self,L):
        out=[]
        for i in range(len(L)):
            if L[i]!=0: # un individu n'a pas de descendant
                out.append(i)
        return(out)
        
    def create_Mating(self,parpop):
        recomb_loci = range(sum(list(parpop.numLoci()))-1) # for the computation of recombination rate per locus we don't take into account the last locus (because it has already been taken into account by the computation of previous loci)
        N_ind = parpop.popSize()
        mating = sim.HeteroMating([
            
            sim.SelfMating(
            ops = [sim.ParentsTagger(),sim.PyTagger(self.autoTagger),sim.Recombinator(
                rates=self.calcRecombination(parpop,recomb_loci),
                loci=recomb_loci)],
            weight = self.selfing,
            selectionField = self.selection),
            
            sim.HermaphroditicMating(
            allowSelfing = False,
            ops = [sim.ParentsTagger(),sim.PyTagger(self.alloTagger),sim.Recombinator(
                rates=self.calcRecombination(parpop,recomb_loci),
                loci=recomb_loci)],
            weight = 1-self.selfing,
            selectionField = self.selection)],
            
            subPopSize = self.newSize)
        return(mating)
    
    def calcRecombination(self,parpop,list_loci):
        """ This function calculates recombination rates for all loci positions """
        rateRecomb = []
        for locus in list_loci:
            chrom1 = parpop.chromLocusPair(locus)[0]
            chrom2 = parpop.chromLocusPair(locus+1)[0]
            if chrom1 == chrom2 :
                delta = parpop.locusPos(locus+1)-parpop.locusPos(locus)
                if self.opt == "Haldane":
                    rate = 0.5*(1- math.exp(-2*delta))
                if self.opt == "Kosambi":
                    rate = np.tanh(2*delta)
                rateRecomb.append(rate)
            else :
                rateRecomb.append(0.5)
        return(rateRecomb)




    
        
    
