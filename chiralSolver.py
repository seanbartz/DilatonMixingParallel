#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 13:38:10 2022
This program moves the solving of the chiral field to its own function.

The main portion of the program loops over T and calculates the values of sigma for each temperature in parallel.

This is NOT faster than looping over sigma and parallelizing over T

@author: seanbartz
"""
import numpy as np
from scipy.integrate import odeint
from solveTmu import blackness

from timebudget import timebudget
import time


from threeflavorALLvalues import chiral
from threeflavorALLvalues import get_all_sigmas_parallel

from multiprocessing import Pool
import os


# @timebudget
def chiralSolver(args):
    zh,q,ml,a0,lambda1,sl=args
    
     
    mu_g=440


    """
    limits of spatial variable z/zh. Should be close to 0 and 1, but 
    cannot go all the way to 0 or 1 because functions diverge there
    """
    ui = 0.01
    uf = 0.999
    "Create the spatial variable mesh"
    umesh=100
    u=np.linspace(ui,uf,umesh)
    
    Q=q*zh**3

    
    
 
    
    "This is a constant that goes into the boundary conditions"
    zeta=np.sqrt(3)/(2*np.pi)
    
    
    "Matching Fang paper"
    v4=4.2
    v3= -22.6/(6*np.sqrt(2))
    
    "need the dilaton for mixing term in test function"
    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)
        
    params=v3,v4,lambda1,mu_g,a0,zh,q
    "blackness function and its derivative, Reissner-Nordstrom metric"
    "This version is for finite temp, finite chemical potential"
    f = 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp = -4*(1+Q**2)*u**3 + 6*Q**2*u**5
    
    "stepsize for search over sigma"
    "Note: search should be done over cube root of sigma, here called sl"
    #tic = time.perf_counter()

    
    
    s2=-3*(ml*zeta)**2*v3
    s3=-9*(zeta*ml)**3*v3**2 + 2*(zeta*ml)**3*v4 + ml*zeta*mu_g**2 - 1/2*ml*zeta*lambda1*mu_g**2
    
    
    # sigmaRange=np.linspace(minsigma,maxsigma,round(maxsigma/deltasig))
    # testIR=np.zeros(len(range (minsigma,maxsigma,deltasig)))

    # for sl in range (minsigma,maxsigma,deltasig):
    "values for chiral field and derivative at UV boundary"
    sigmal = sl**3
    UVbound = [ml*zeta*zh*ui + sigmal/zeta*(zh*ui)**3+s2*(zh*ui)**2+s3*(zh*ui)**3*np.log(zh*ui), 
               ml*zeta*zh + 3*sigmal/zeta*zh**3*ui**2 + 2*s2*zh**2*ui + s3* ui**2*zh**3*(1+3*np.log(zh*ui))]
       
    "solve for the chiral field"
    chiFields=odeint(chiral,UVbound,u,args=(params,))
    
    "test function defined to find when the chiral field doesn't diverge"
    "When test function is zero at uf, the chiral field doesn't diverge"
    test = ((-u**2*fp)/f)*chiFields[:,1]-1/f*(3*chiFields[:,0]+lambda1*phi*chiFields[:,0]-3*v3*chiFields[:,0]**2-4*v4*chiFields[:,0]**3)
    testIR = test[umesh-1]#value of test function at uf
    return testIR

if __name__ == '__main__':
    tmin=102
    tmax=102.5
    numtemp=10

    temps=np.linspace(tmin,tmax,numtemp)
    ml=24
    
    #chemical potential
    mu=500
    
    lambda1= 7.438#parameter for mixing between dilaton and chiral field

    
    a0=0.

    
    minsigma=0
    maxsigma=300
    deltasig=1
    
    sigmaRange=np.linspace(minsigma,maxsigma,round(maxsigma/deltasig))
    testIR=np.zeros([numtemp,len(range (minsigma,maxsigma,deltasig))])
    
    
    processes_count=os.cpu_count()    
    processes_pool2 = Pool(processes_count)
    
    
    start_time=time.perf_counter()
    for i in range(numtemp):
        "solve for horizon and charge"
        args=temps[i],mu,ml,a0,lambda1

        # zh,q=blackness(temps[i],mu)
        argsArray=np.vstack([np.outer(args,np.ones(len(sigmaRange))),sigmaRange]).T
    
        testIR[i,:]=processes_pool2.map(chiralSolver,argsArray)
        
    print("process takes", time.perf_counter()-start_time, " seconds")
    processes_pool2.close()
    
    # print(chiralSolver([zh,q,ml,a0,lambda1,240]))
    