#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

EDITED on March 10 2023:
    Allows stepsize over sigma to be non-integer, for finer search

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
from scipy.integrate import odeint
from solveTmu import blackness

from timebudget import timebudget

import matplotlib.pyplot as plt
from multiprocessing import Pool
import os


# import time




# start_time=time.perf_counter()

def chiral(y,u,params):
    chi,chip=y
    v3,v4,lambda1,mu_g,a0,zh,q=params
    
    Q=q*zh**3
    
    
    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)
    phip = 2*u*(zh*mu_g)**2+a0*(4*u**6*(zh*mu_g)**7/(1+(u*zh*mu_g)**4)**2-3*u**2*(zh*mu_g)**3/(1+(u*zh*mu_g)**4))

    f= 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp= -4*(1+Q**2)*u**3 + 6*Q**2*u**5
    "EOM for chiral field"
    derivs=[chip,
            (3/u-fp/f+phip)*chip - (3*chi+lambda1*phi*chi-3*v3*chi**2-4*v4*chi**3)/(u**2*f)]
            #((3+u**4)/(u-u**5) +phip)*chip - (-3*chi+4*v4*chi**3)/(u**2-u**6) ]
            
    return derivs
# @timebudget
def allSigmas(args):#,mu,ml,minsigma,maxsigma,a0,lambda1):
    "Unpack the input"
    T,mu,ml,minsigma,maxsigma,a0,lambda1=args

    minsigma=int(minsigma)
    maxsigma=int(maxsigma)
    "stepsize for search over sigma"
    "Note: search should be done over cube root of sigma, here called sl"
    deltasig = 0.5
    
    # create an array of sigma values from minsigma to maxsigma, incrementing by deltasig
    sigmavalues = np.arange(minsigma,maxsigma,deltasig)
    
    mu_g=440

    "solve for horizon and charge"
    zh,q=blackness(T,mu)
    Q=q*zh**3
    """
    limits of spatial variable z/zh. Should be close to 0 and 1, but 
    cannot go all the way to 0 or 1 because functions diverge there
    """
    ui = 1e-2
    uf = 1-ui
    "Create the spatial variable mesh"
    umesh=100
    u=np.linspace(ui,uf,umesh)
    

    
    
 
    
    "This is a constant that goes into the boundary conditions"
    zeta=np.sqrt(3)/(2*np.pi)
    
    "For the scalar potential in the action"
    "see papers by Bartz, Jacobson"
    #v3= -3 #only needed for 2+1 flavor
    # v4 = 8
    # v3 = -3
    
    "Matching Fang paper"
    v4=4.2
    v3= -22.6/(6*np.sqrt(2))
    
    "need the dilaton for mixing term in test function"
    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)
        
    #sigmal=260**3
    params=v3,v4,lambda1,mu_g,a0,zh,q
    "blackness function and its derivative, Reissner-Nordstrom metric"
    "This version is for finite temp, finite chemical potential"
    f = 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp = -4*(1+Q**2)*u**3 + 6*Q**2*u**5
    

    #tic = time.perf_counter()

    truesigma = 0
    "This version steps over all values to find multiple solutions at some temps"
    
    "initial values for comparing test function"
    oldtest=0
    j=0
    truesigma=np.zeros(3)
    
    
    s2=-3*(ml*zeta)**2*v3
    s3=-9*(zeta*ml)**3*v3**2 + 2*(zeta*ml)**3*v4 + ml*zeta*mu_g**2 - 1/2*ml*zeta*lambda1*mu_g**2

    #use this line if deltasig is an integer
    # for sl in range (minsigma,maxsigma,deltasig):
    
    #use these next two lines if deltasig is not an integer
    for i in range(len(sigmavalues)):
        sl=sigmavalues[i]
    
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
        
        "when test function crosses zero, it will go from + to -, or vice versa"
        "This is checked by multiplying by value from previous value of sigma"
        if oldtest*testIR<0: #and chiFields[umesh-1,0]>0:
           
            truesigma[j]=sl #save this value
            j=j+1 #if there are other sigma values, they will be stored also
            #print(truesigma)
        if j>2:
            break
            
        oldtest=testIR

    
    return truesigma

@timebudget
def get_all_sigmas(operation, input):
    "This function executes a loop to calculate all sigma values for all values of the temps array"
    truesigma=np.zeros([len(input),3])

    for i in range(0,len(input)):
        truesigma[i,:]=operation(input[i])#,100,24,0,300,0,7.438)
    return truesigma

@timebudget
def get_all_sigmas_parallel(operation,input,pool):
    truesigma=np.zeros([len(input),3])

    truesigma=pool.map(operation, input)
    
    return truesigma

    

if __name__ == '__main__':

        
    tmin=50
    tmax=100

    numtemp=100
    
    temps=np.linspace(tmin,tmax,numtemp)
    
    #light quark mass
    ml=0*np.ones(numtemp)
    
    #chemical potential
    mu=0*np.ones(numtemp)
    
    lambda1=5.5*np.ones(numtemp) #parameter for mixing between dilaton and chiral field
    
    minsigma=0*np.ones(numtemp)
    maxsigma=200*np.ones(numtemp)
    
    a0=0.*np.ones(numtemp)
    
    
    tempsArgs=np.array([temps,mu,ml,minsigma,maxsigma,a0,lambda1]).T


    #need up to 3 sigma values per temperature
    # truesigma=np.zeros([numtemp,3])
    
    "This calls the old version, which loops over all temps. Only un-comment for speed comparisons"
    # truesigma=get_all_sigmas(allSigmas,tempsArgs)
    
    "Create a pool that uses all available cpus"
    processes_count=os.cpu_count()    
    processes_pool = Pool(processes_count)
    
    truesigma=get_all_sigmas_parallel(allSigmas,tempsArgs,processes_pool)
    truesigma=np.array(truesigma)**3/1e9
    processes_pool.close()
    
    # find the indices of the non-zero values of truesigma
    # this is necessary because the function returns 0 for sigma values that don't exist
    nonzero1=np.nonzero(truesigma[:,0])
    truesigma1=truesigma[:,0][nonzero1]
    temps1=temps[nonzero1]

    nonzero2=np.nonzero(truesigma[:,1])
    truesigma2=truesigma[:,1][nonzero2]
    temps2=temps[nonzero2]

    nonzero3=np.nonzero(truesigma[:,2])
    truesigma3=truesigma[:,2][nonzero3]
    temps3=temps[nonzero3]

    # scatter plot the non-zero values    
    plt.scatter(temps1,truesigma1)
    # plt.scatter(temps2,truesigma2)
    # plt.scatter(temps3,truesigma3)
    #plt.ylim([min(truesigma1)-5,max(truesigma[:,0])+5])
    plt.xlabel('Temperature (MeV)')
    plt.ylabel(r'$\sigma^{1/3}$ (MeV)')
    plt.title(r'$m_q=%i$ MeV, $\mu=%i$ MeV, $\lambda_1=$ %f' %(ml[0],mu[0],lambda1[0]))
    plt.show()

    if max(truesigma[:,1])==0:
        print("Crossover or 2nd order")
        #find the temp value where the gradient of truesigma[:,0] is most negative
        #this is the pseudo-critical temperature
        print("Pseudo-Critical temperature is between", temps[np.argmin(np.gradient(truesigma[:,0]))-1], temps[np.argmin(np.gradient(truesigma[:,0]))] )
        #these temperature values are the new bounds for the next iteration
        tmin=temps[np.argmin(np.gradient(truesigma[:,0]))-1]
        tmax=temps[np.argmin(np.gradient(truesigma[:,0]))]

        #these values of sigma are the new bounds for the next iteration
        maxsigma=truesigma[np.argmin(np.gradient(truesigma[:,0]))-1,0]
        minsigma=truesigma[np.argmin(np.gradient(truesigma[:,0])),0]

        #print the sigma values for the new bounds
        print("Sigma bounds for the next search are ", minsigma, maxsigma)

        #plot just truesigma1
        plt.plot(temps1,truesigma1,linewidth=3)
        plt.xlabel('Temperature (MeV)')
        plt.ylabel(r'$\sigma^{1/3}$ (MeV)')
        plt.title(r'$m_q=%i$ MeV, $\mu=%i$ MeV, $\lambda_1=$ %f' %(ml[0],mu[0],lambda1[0]))
        plt.show()
    else:
        print("First order")  
        #crtical temperature is where truesigma has multiple solutions  
        print("Critical temperature is ", temps[np.argmax(truesigma[:,1])] )

        #reverse the order of the arrays where the plot goes "backward"
        temps2= temps2[::-1]
        truesigma2= truesigma2[::-1]
        #find the index where the gradient of truesgima1 is most negative
        splitIndex1=np.argmin(np.gradient(truesigma1))

        #split truesigma1 and temps1 into two arrays at the index where the gradient is most negative
        truesigma1a=truesigma1[:splitIndex1]
        truesigma1b=truesigma1[splitIndex1:]
        temps1a=temps1[:splitIndex1]
        temps1b=temps1[splitIndex1:]

        #join the arrays in the following order truesigma1a, truesigma3, truesigma2, truesigma1b
        truesigma=np.concatenate((truesigma1a,truesigma3,truesigma2,truesigma1b))
        temps=np.concatenate((temps1a,temps3,temps2,temps1b))

        "Uncomment these lines to make a nice-looking graph"
        plt.plot(temps,truesigma,linewidth=3)
        plt.xlabel('Temperature (MeV)')
        plt.ylabel(r'$\sigma^{1/3}$ (MeV)')
        plt.title(r'$m_q=%i$ MeV, $\mu=%i$ MeV, $\lambda_1=%.3f$' %(ml[0],mu[0],float(str(lambda1[0]).rstrip('0').rstrip('.'))))        
        # save the plot 
        plt.savefig('plots/sigmaVsT_lambda_7438_mq_24_mu_400.png',dpi=500)
        
        plt.show()
    # end_time=time.perf_counter()
    # print("Time elapsed = ", end_time-start_time )

    #export temps and truesigma to a the same csv file for plotting in another program
    # np.savetxt('sigmaVsT_lambda_6_mq_15_mu_500.csv',np.column_stack((temps,truesigma)),delimiter=',',comments='')
        