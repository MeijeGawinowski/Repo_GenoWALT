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
    
    def __init__(self,N):
        self.generations = 1
        self.replicat = 1
        self.selfing = 0.5 # proportion of autogamy in the population
        self.selection = 'fitness' # selection trait for reproduction
        self.opt = "Haldane" # computation option for recombination rates
        self.newSize = N
        self.attr = "drawing"

            
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
        
        Des = list(parpop.indInfo('fitness'))
                
        N_ind = parpop.popSize()
        ind_idx=list(range(N_ind))
        
        DesTot=sum(Des)
        prop=[]
        statN_des=[]
        if self.attr == "mult":
            for i in range(len(Des)):
                p = Des[i]/DesTot
                prop.append(p)
                n = round(self.newSize*p)
                statN_des.append(n)
        if self.attr == "drawing" :
            for i in range(len(Des)):
                p = Des[i]/DesTot
                prop.append(p)
            par = list(np.random.choice(list(range(N_ind)),size=self.newSize, replace=True, p=prop))
            statN_des = [par.count(i) for i in range(N_ind)]
        # print("Nb of offspring per parent : ",statN_des)
        dynN_des=cp.deepcopy(statN_des) # updated during reproduction process    
        newpop=cp.deepcopy(parpop)
        newpop.removeIndividuals(indexes=range(N_ind))
        newpop.removeInfoFields(list(parpop.infoFields()))
        newpop.addInfoFields(['mode','mother_idx','father_idx'])
        newpop.setSubPopName('test',0)
            
        for i in range(N_ind):
            # for the individual i
            # print("Progeny of individual N°",i)
            scheme=list(np.random.choice([0,1],size=dynN_des[i], replace=True,p=[self.selfing,1-self.selfing]))
            # print(scheme)
            N_auto=scheme.count(0)
            N_allo=scheme.count(1)
            # print("allof=",N_allo," and autof=",N_auto)
                
            # production of offsprings by allofecondation
            # print("ALLOF")
            for j in range(N_allo): # for the j offsprings produced by allofecondation by individual i
                list_par = list(range(N_ind))
                copy_prop = cp.deepcopy(prop)
                del list_par[i] # delete the parent i itself (no self-fertilization)
                del copy_prop[i]
                Sprop = sum(copy_prop)
                copy_prop = [copy_prop[i]/Sprop for i in range(len(copy_prop))]
                lover=int(np.random.choice(list_par,size=1,replace=True,p=copy_prop))
                # print("Partner nb",j," for ",i," is ind nb",lover)
                    
                    # dynN_des[lover]=dynN_des[lover]-1
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
            
            
            # print("AUTOF")
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




    
        
    
