#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 17:29:19 2022
Parallelize using pool
@author: seanbartz
"""

import numpy as np
# from scipy.integrate import odeint

# from solveTmu import blackness
from threeflavorRefinedmu import sigmasearch

from timebudget import timebudget
from multiprocessing import Pool
import os



"light quark mass"
quarkmass=24
"quark chemical potential"
chemPotential=100

lambda1=7.438
a0=0

@timebudget
def get_sigmas(operation, input):
    "This function executes a loop to calculate sigma for all values of the temps array"
    sigmaArray=np.zeros(len(input))
    for i in range(0,len(input)):
        sigmaArray[i]=operation(input[i])#,100,24,7.438,0)
    return sigmaArray
        
@timebudget
def parallel_get_sigmas(operation,input,pool):
    "This function uses pools to calculate sigm for all values of temps in parallel"
    # args=[100,24,7.438,0]
    sigmaArray=np.zeros(len(input))

    sigmaArray=pool.map(operation, input)
    return sigmaArray
    
            
        
if __name__ == '__main__':
    mintemp=100
    maxtemp=150

    perMev=1 #number of data points per MeV of temperature (must be integer)
    numtemps=(maxtemp-mintemp)*perMev+1
    temps=np.linspace(mintemp,maxtemp,numtemps)
    sigmas=get_sigmas(sigmasearch,temps)
    
    "Create a pool that uses all available cpus"
    processes_count=os.cpu_count()    
    processes_pool = Pool(processes_count)
    
    sigmas2=parallel_get_sigmas(sigmasearch,temps,processes_pool)
    processes_pool.close()

    



