#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EDITED ON June 9 2022:
    Added mixing term between dilaton and chiral field
Created on Tue Feb 9  2021
Uses a refining search process to speed up search
Finds the correct value of sigma for a given temperature
Flavor-symmetric (three flavor)
Finite temperature, finite chemical potential 

Set up to be called as outside function of temperature, chem potential, and light quark mass
@author: seanbartz
"""
import numpy as np
from scipy.integrate import odeint
from solveTmu import blackness
from timebudget import timebudget


#import time







# "temperature in MeV"
# Temp=150
# #light quark mass
# ml=30




def chiral(y,u,params):
    chi,chip=y
    v3,v4,lambda1,mu_g,a0,zh,q=params
    
    Q=q*zh**3

    
    # phi = -(mu1*zh*u)**2 + (mu1**2+mu0**2)*(zh*u)**2*(1 - np.exp(-(mu2*zh*u)**2))
    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)
    phip = 2*u*(zh*mu_g)**2+a0*(4*u**6*(zh*mu_g)**7/(1+(u*zh*mu_g)**4)**2-3*u**2*(zh*mu_g)**3/(1+(u*zh*mu_g)**4))
    
    "derivative of the dilaton, using exp parameterization"
    # phip= 2*u*zh**2*(mu0**2+np.exp(-(mu2*zh*u)**2)*(mu0**2+mu1**2)*((u*zh*mu2)**2-1) )
    "blackness function and its derivative, Reissner-Nordstrom metric"
    "This version is for finite temp, finite chemical potential"
    f= 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp= -4*(1+Q**2)*u**3 + 6*Q**2*u**5
    "EOM for chiral field"
    derivs=[chip,
            (3/u-fp/f+phip)*chip - (3*chi+lambda1*phi*chi-3*v3*chi**2-4*v4*chi**3)/(u**2*f)]
            #((3+u**4)/(u-u**5) +phip)*chip - (-3*chi+4*v4*chi**3)/(u**2-u**6) ]
            
    return derivs

# @timebudget
def sigmasearch(T):#,mu,ml,lambda1,a0):
    mu=100
    ml=24
    lambda1=7.438
    a0=0
    "solve for horizon and charge"
    zh,q=blackness(T,mu)
    Q=q*zh**3
    """
    limits of spatial variable z/zh. Should be close to 0 and 1, but 
    cannot go all the way to 0 or 1 because functions diverge there
    """
    ui = 0.01
    uf = 0.999
    "Create the spatial variable mesh"
    umesh=100
    u=np.linspace(ui,uf,umesh)
    
    
    #parameters for dilaton. See papers
    mu_g=440

    
    
 
    
    "This is a constant that goes into the boundary conditions"
    zeta=np.sqrt(3)/(2*np.pi)
    
    "For the scalar potential in the action"
    "see papers by Bartz, Jacobson"
    #v3= -3 #only needed for 2+1 flavor
    #v4 = 8
    #v3 = -3
    "Matching Fang paper"
    v4=4.2
    v3= -22.6/(6*np.sqrt(2))
        
    #sigmal=260**3
    params=v3,v4,lambda1,mu_g,a0,zh,q
    "blackness function and its derivative, Reissner-Nordstrom metric"
    "This version is for finite temp, finite chemical potential"
    f = 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp = -4*(1+Q**2)*u**3 + 6*Q**2*u**5
    

    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)

    "stepsize for search over sigma"
    "Note: search should be done over cube root of sigma, here called sl"
    deltasig = 100
    #tic = time.perf_counter()
    minsigma = 0
    maxsigma = 500
    "initialize variable for correct sigma value"
    truesigma = 0#[]
    "This version uses a refining method to search"
    "Runs an order of magnitude faster"
    s2=-3*(ml*zeta)**2*v3
    s3=-9*(zeta*ml)**3*v3**2 + 2*(zeta*ml)**3*v4 + ml*zeta*mu_g**2 - 1/2*ml*zeta*lambda1*mu_g**2
    while deltasig > 0.1:
    
        "initial values for comparing test function"
        oldtest=0
        #print(deltasig)


        for sl in range (minsigma,maxsigma,deltasig):
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
            
            "Breaks the loop if the chiral field is blowing up"
            "This is only here to save time, so comment out if it causes problems"
            if np.amax(abs(chiFields[:,0]))>1e5:
                deltasig=0.01
                break
            
            "when test function crosses zero, it will go from + to -, or vice versa"
            "This is checked by multiplying by value from previous value of sigma"
            if oldtest*testIR<0: #and chiFields[umesh-1,0]>0:
                #print(oldtest*testIR)
                #print(sl)
                #print(chiFields[umesh-1,0])
                truesigma=sl #save this value
                "new range is +/- deltasig"
                "Need + and - in case true value is right on a multiple of deltasig"
                maxsigma=sl+deltasig
                minsigma=sl-deltasig
                deltasig=int(deltasig*.1)
                break
            oldtest=testIR
        #print(sl)
        "This protects the program from getting hung up"
        "if it reaches the top of the search range, it refines the search"
        "It may not find anything there, but it will eventually terminate"
        if sl >= maxsigma-deltasig:
            deltasig=int(deltasig*.1)
    #toc = time.perf_counter()
    #print(f"Found the value in {toc - tic:0.4f} seconds")
    "solve for the chiral field with the correct physical value of sigma"
    
    return truesigma#,np.amax(abs(chiFields[:,0]))
'call the function with arguments temperature, chemical potential,  quark mass, and dilaton-chiral mixing parameter'
if __name__ == '__main__':
    # import matplotlib.pyplot as plt

    (truesigma,)=sigmasearch(100)#,100,24,7.438,0)
    print("sigma value is ", truesigma)

    # sigmal=sl**3
    # UVbound =[ml*zeta*zh*ui + sigmal/zeta*(zh*ui)**3, ml*zeta*zh + 3*sigmal/zeta*zh**3*ui**2]
    # chiFields=odeint(chiral,UVbound,u,args=(params,))    
    # test=((-u**2*fp)/f)*chiFields[:,1]-1/f*(3*chiFields[:,0]-4*v4*chiFields[:,0]**3)

    # "plot the chiral field and test function for the correct sigma value"
    # fig, ax1=plt.subplots()
    # plt.plot(u,chiFields[:,0])
    # plt.plot(u,test)
    # plt.plot(u,ml*zeta*zh*u + sigmal/zeta*(zh*u)**3)
    # plt.show
