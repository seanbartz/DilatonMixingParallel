#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 15:48:44 2022

Figured out how to pass additional parameters while using Pools.

-------
Branched from earlier file (threeflavorALLvalues.py) with the following notes:

EDITED on November 15 2022:
    Parallel processing with Pools for speedup by factor of ~3 on an 8-CPU machine

EDITED ON June 9 2022:
    Added mixing term between dilaton and chiral field
    
Created on Tue March 16  2021
For a given quark mass and chemical potential, 
solves for all sigma values for a range of temperatures.
If there are multiple values, then the transition is 1st order.

@author: seanbartz
"""

import numpy as np
from solveTmu import blackness
from threeflavorALLvalues import chiral

from scipy.integrate import odeint
from timebudget import timebudget

import matplotlib.pyplot as plt
from multiprocessing import Pool
import os

from threeflavorALLvalues import get_all_sigmas_parallel

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

# @timebudget
def allSigmas(args):
    "Unpack the input"
    T,mu,ml,minsigma,maxsigma,a0,lambda1=args

    "solve for horizon and charge"
    zh,q=blackness(T,mu)
    

    minsigma=int(minsigma)
    maxsigma=int(maxsigma)
    
    "stepsize for search over sigma"
    "Note: search should be done over cube root of sigma, here called sl"
    deltasig = 1
    #tic = time.perf_counter()

    "This version steps over all values to find multiple solutions at some temps"
    

    
    sigmaRange=np.linspace(minsigma,maxsigma,round(maxsigma/deltasig))
    testIR=np.zeros(len(range (minsigma,maxsigma,deltasig)))
    
    '''
    This was intended to parallelize the solving the chiral field at various values of sigma,
    however, this is not allowed, as 'daemonic processes are not allowed to have children'
    
    # processes_count=os.cpu_count()    
    # processes_pool2 = Pool(processes_count)
    
    # argsArray=np.vstack([np.outer(args,np.ones(len(sigmaRange))),sigmaRange]).T
    
    # testIR=chiralSolver(argsArray,processes_pool2)
    # processes_pool2.close()
    '''
    
    
    i=0
    for sl in range (minsigma,maxsigma,deltasig):
        "values for chiral field and derivative at UV boundary"
        testIR[i] = chiralSolver([zh,q,ml,a0,lambda1,sl])#value of test function at uf
        i=i+1

    signChange = np.where(testIR[:-1] * testIR[1:] < 0 )[0] +1
    truesigma=sigmaRange[signChange]
    truesigma=truesigma[0:3]
    truesigma=np.pad(truesigma,(0,3-len(truesigma)),'constant')
    
    return truesigma
if __name__ == '__main__':
    
    tmin=102
    tmax=102.5
    numtemp=10

    temps=np.linspace(tmin,tmax,numtemp)
    
        #light quark mass
    ml=24*np.ones(numtemp)
    
    #chemical potential
    mu=500*np.ones(numtemp)
    
    lambda1= 7.438*np.ones(numtemp) #parameter for mixing between dilaton and chiral field
    
    minsigma=00*np.ones(numtemp)
    maxsigma=300*np.ones(numtemp)
    
    a0=0.*np.ones(numtemp)
    
    
    tempsArgs=np.array([temps,mu,ml,minsigma,maxsigma,a0,lambda1]).T
    "Create a pool that uses all available cpus"
    processes_count=os.cpu_count()    
    processes_pool = Pool(processes_count)
    
    truesigma=get_all_sigmas_parallel(allSigmas,tempsArgs,processes_pool)
    truesigma=np.array(truesigma)
    processes_pool.close()
    
        
    plt.scatter(temps,truesigma[:,0])
    plt.scatter(temps,truesigma[:,1])
    plt.scatter(temps,truesigma[:,2])
    plt.ylim([min(truesigma[:,0])-5,max(truesigma[:,0])+5])
    plt.xlabel('Temperature (MeV)')
    plt.ylabel(r'$\sigma^{1/3}$ (MeV)')
    plt.title(r'$m_q=%i$ MeV, $\mu=%i$ MeV, $\lambda_1=$ %f' %(ml[0],mu[0],lambda1[0]))
    
    #save the figure with the parameters in the name, in a subfolder called 'plots'
    plt.savefig('plots/sigma_vs_T_mq%i_mu%i_lam%f.eps' %(ml[0],mu[0],lambda1[0]))

    #save the data to a npy file with the parameters in the name, in a subfolder called 'data'
    np.save('data/sigma_vs_T_mq%i_mu%i_lam%f.npy' %(ml[0],mu[0],lambda1[0]),truesigma)

    if max(truesigma[:,1])==0:
        print("Crossover or 2nd order")
    else:
        print("First order")    
        print("Critical temperature is ", temps[np.argmax(truesigma[:,1])] )