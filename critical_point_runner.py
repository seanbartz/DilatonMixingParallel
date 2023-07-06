import sys
import critical_point

lambda1 = float(sys.argv[1])
a0 = float(sys.argv[2])
ml = float(sys.argv[3])
tmin = float(sys.argv[4])
tmax = float(sys.argv[5])
numtemp = int(sys.argv[6])
minsigma = float(sys.argv[7])
maxsigma = float(sys.argv[8])
mu_initial = float(sys.argv[9])
delta_mu = float(sys.argv[10])
mu_precision = int(sys.argv[11])

critical_point.critical_point_refined(lambda1,a0,ml,tmin,tmax,numtemp,minsigma,maxsigma,mu_initial,delta_mu,mu_precision)
