import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

from scipy.optimize import brentq
from scipy.optimize import newton

from threeflavorALLvalues import chiral
from solveTmu import blackness


#calculates the value of the test function at the horizon for a given sigma and other parameters
#returns the value of the test function at the horizon
def testHorizon(sl,ml,T, mu,mu_g,a0,lambda1):


    #sl,ml,T,mu_g,a0,lambda1=args

    "solve for horizon and charge"
    zh,q=blackness(T,mu)
    Q=q*zh**3    # zi=0.001/(np.pi*140) #try fixing zi to be a constant based on T=140
    ui =0.001# zi/zh #0.001
    uf = .9
    "Create the spatial variable mesh"
    umesh=100
    u=np.linspace(ui,uf,umesh)

    
    "This is a constant that goes into the boundary conditions"
    zeta=np.sqrt(3)/(2*np.pi)
    
        
    params=lambda1,mu_g,a0,zh
    "blackness function and its derivative, Schwarzschild metric"
    "This version is for finite temp, zero chemical potential"
    f= 1 - (1+Q**2)*u**4 + Q**2*u**6
    fp= -4*(1+Q**2)*u**3 + 6*Q**2*u**5


    "Matching Fang paper"
    v4=4.2
    v3= -22.6/(6*np.sqrt(2))
    
    "need the dilaton for mixing term in test function"
    "Ballon-Bayona version"
    phi = (mu_g*zh*u)**2-a0*(mu_g*zh*u)**3/(1+(mu_g*zh*u)**4)

    sigmal = sl**3
    params=v3,v4,lambda1,mu_g,a0,zh,q

    s2=-3*(ml*zeta)**2*v3
    s3=-9*(zeta*ml)**3*v3**2 + 2*(zeta*ml)**3*v4 + ml*zeta*mu_g**2 - 1/2*ml*zeta*lambda1*mu_g**2


    UVbound = [ml*zeta*zh*ui + sigmal/zeta*(zh*ui)**3+s2*(zh*ui)**2+s3*(zh*ui)**3*np.log(zh*ui), 
                   ml*zeta*zh + 3*sigmal/zeta*zh**3*ui**2 + 2*s2*zh**2*ui + s3* ui**2*zh**3*(1+3*np.log(zh*ui))]
    "solve for the chiral field"
    chiFields=odeint(chiral,UVbound,u,args=(params,))
    #chi=chiFields[:,0]
    "solve for the metric"

    "test function defined to find when the chiral field doesn't diverge"
    "When test function is zero at uf, the chiral field doesn't diverge"
    test = ((-u**2*fp)/f)*chiFields[:,1]-1/f*(3*chiFields[:,0]-3/2*lambda1*chiFields[:,0]**3)
    testIR = test[umesh-1]#value of test function at uf

    return testIR

# Find one to three roots of testHorizon in the interval [minsigma, maxsigma]
# using the Brent method.
# The function testHorizon is called with the arguments (sigma, ml, T, mu_g, a0, lambda1).
# The roots are returned as a list of tuples (sigma, testHorizon(sigma)).
# The number of roots found is returned as an integer.
def findRootsBrent(minsigma, maxsigma, ml, T, mu, mu_g, a0, lambda1):
    #Note: the brent method only finds one root at a time, 
    # and only finds roots if the function values at the endpoints have opposite signs.
    # So we have to call the function three times to find three roots, and guess the intervals containing the roots
    roots = []

    #find the first  root
    try:
        sigma = newton(testHorizon,  maxsigma, args=(ml, T, mu, mu_g, a0, lambda1))
        roots.append(sigma)
    except ValueError:
        pass

    #find the second root if a first root was found
    # in the interval [minsigma/2,sigma]


    try:
        sigma = brentq(testHorizon, minsigma/2, sigma, (ml, T, mu, mu_g, a0, lambda1))
        roots.append(sigma)
    except ValueError:
        pass

    #find the third root if a second root was found
    # in the interval [minsigma,minsigma/2]

    try:
        sigma = brentq(testHorizon, minsigma, minsigma/2, (ml, T, mu, mu_g, a0, lambda1))
        roots.append(sigma)
    except ValueError:
        pass
    numRoots=len(roots)

    #if len roots <3, pad the list with zeros
    if numRoots < 3:
        roots = roots + [0]*(3-numRoots)


    return numRoots, roots

#create the main function
if __name__ == '__main__':
    minsigma = 0
    maxsigma = 300

    mu_g = 440

    a0 = 0
    lambda1 = 7.835
    T = 150
    mu= 10
    ml = 10

    #find the roots
    nroots, roots = findRootsBrent(minsigma, maxsigma, ml, T, mu,mu_g, a0, lambda1)
    
    
    
    print(nroots,roots)