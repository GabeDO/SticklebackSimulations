# -*- coding: utf-8 -*-
"""
Created on Thu Mar 24 12:13:03 2022

@author: Gabe
"""


from simuPOP import *
from simuPOP.utils import saveCSV
from simuPOP.demography import *
import simuPOP as sim
import time
import csv
from simuPOP.utils import importPopulation, export



print("Reading Allele Frequencies")


#Gen0 Allele Freqs.csv
#Import alll the Allele Fequencies for each loci
AlleleFreqData = 'INPUT DATA NAME.csv'

with open(AlleleFreqData, newline = '') as f:
    reader = csv.reader(f)
    AlleleFrequencyList = list(reader)
print("Finished Reading Allele Frequencies")


#set somne parameters to be used later
LocNumber = len(AlleleFrequencyList)
print("number of loci: ", LocNumber)
start = time.time()


print("Setting Up Header For CSV File")
#Set up the header for the output file
with open('9 gens - gen every 1 year1 (Model2).csv', 'a', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     header = ["Replicate", "Generation", "popSize"]
     for locus in AlleleFrequencyList:
         header.append(str(locus[0])+"-"+str(locus[1]))
     wr.writerow(header)
print("Header Is Set Up") 



# a little function to keep a track of how long everything is going
def timecount(pop,param):    
    end = time.time()   
    time_elapsed = end - param
    print(time_elapsed)
    return True

#a function to export data - applied each generations
def exportFunc(pop, param):
    
    #tell the code to get allele freq data
    stat(pop, alleleFreq=ALL_AVAIL, popSize = True, subPops = ALL_AVAIL)
    Value = []
    #make a valuye with the repliocate (param) and the generation (pop.dvars().gen+1). It's +1 because python counts from 0
    Info = [param,pop.dvars().gen+1,pop.dvars().subPopSize[0]]
    
    # now run through each loci and 
    for i in pop.dvars().alleleFreq:
        #make a list of all the allele freq per loci
        Value.append(pop.dvars().alleleFreq[i][0])
        
        #if you want expected het, use this:
        #Value.append(1 - sum([x*x for x in pop.dvars().alleleFreq[i].values()]))
    
    #write all this to a new line in the output data
    with open('9 gens - gen every 1 year1 (Model2).csv', 'a', newline='') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(Info + Value)
        
    return True

#1st model with a generation time of 2 years (so 5 generations total in the sim)
model1 = EventBasedModel(
    N0 = (1933, 'Gen 0 '),
    T = 5,
    events=[
        ResizeEvent(at = 1, sizes = 1656),
        ResizeEvent(at = 2, sizes = 1410),
        ResizeEvent(at = 3, sizes = 1956),
        ResizeEvent(at = 4, sizes = 1806),
        ])

#2nd model with a generation time of 1 years (so 9 generations total in the sim), with the pop size changing every 2nd year
model2 = EventBasedModel(
    N0 = (1933, 'Gen 0'),
    T = 9,
    events=[
        ResizeEvent(at = 2, sizes = 1656),
        ResizeEvent(at = 4, sizes = 1410),
        ResizeEvent(at = 6, sizes = 1956),
        ResizeEvent(at = 8, sizes = 1806),
        ])







#set replicates
for rep in range(10):
    
    print("Setting up population for Replicate", str(rep))
    pop = sim.Population(size = model2.init_size, loci=LocNumber)
    stat(pop, alleleFreq=ALL_AVAIL) 

    #here I'm setting the allele freq for each loci for the input population
    lociID = 0
    for i in AlleleFrequencyList:    
        sim.initGenotype(pop, freq=[float(i[2]), 1-float(i[2])], loci=lociID)
        lociID += 1


    print("Population is set up!")
    
    #now we start simulating
    pop.evolve(
        
            #This sim uses random mating
            matingScheme = RandomMating(subPopSize=model2),   
            #and an even sex ratio
            initOps = [
                        InitSex()],
            #here is where we tell it to run our export code each generation, as well as out time keeping code
            postOps = [
                PyOperator(func=timecount, param = start),
                PyOperator(func=exportFunc, param = rep),
            ],
            #we run it for a number of generations specified int he model
            gen=model2.num_gens
            )
    
#and this is just to keep the terminal open after its donem untill I want to exit. 
answer = input("Exit?")