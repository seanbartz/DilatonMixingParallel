import numpy as np
import matplotlib.pyplot as plt
#import odeint
from scipy.integrate import odeint
# import brentq
from scipy.optimize import brentq
from solveTmu import blackness
# import solve_ivp
from scipy.integrate import solve_ivp


def chiral(u,y,params):
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

def chiral_solve_IR(d0,lambda1,T,mu,ui,uf):
    u=np.linspace(ui,uf,1000)
    "Matching Fang paper"
    v4=4.2
    v3= -22.6/(6*np.sqrt(2))

    lambda3=v3
    lambda4=8*v4/3


    zh,q = blackness(T,mu)
    Q=q*zh**3

    d1 = (3 * d0 - 3 * d0**2 * lambda3 - 4 * d0**3 * lambda4 + d0 * zh**2 * lambda1 * mu_g**2) / (2 * (-2 + Q**2))


    "IR boundary condition"
    chi0=d0+d1(1-uf)
    chip0=d1
    y0=[chi0,chip0]

    mu_g=440
    a0=0



    params=v3,v4,lambda1,mu_g,a0,zh,q
    
    "solve the EOM using solve_ivp"
    sol=solve_ivp(chiral,[uf,ui],y0,t_eval=u,args=(params,))

    #sol=odeint(chiral,y0,u,args=(params,))
    #plot
    chi=sol.y[0]
    chip=sol.y[1]
    plt.plot(sol.t,chi)
    plt.show()
    
    return chi,chip,u

if __name__ == "__main__":
    lambda1=7
    T=120
    mu=0
    ui=0.01
    uf=0.999
    d0=0.2
    chi,chip,u=chiral_solve_IR(d0,lambda1,T,mu,ui,uf)

 
