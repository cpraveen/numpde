import os
import numpy as np
from matplotlib import pyplot as plt

# path to executable
exe = '../src/advect'
# check that executable is present
if os.path.isfile(exe)==False:
    print "Could not find ", exe
    exit()

N = [40, 60, 80, 100, 120, 150]
itmax = 50000
itsave = 100
fluxtype = "mda"
limtype  = "weno5"
method = "method6"
cfl = 0.5
testcase = 1

os.system("rm -f err.dat")
for n in N:
    print "n =",n
    # create input file
    f = open('inp.dat','w')
    f.write(str(n)+'\n')
    f.write(str(n)+'\n')
    f.write(str(itmax)+'\n')
    f.write(str(itsave)+'\n')
    f.write(fluxtype+'\n')
    f.write(limtype+'\n')
    f.write(method+'\n')
    f.write(str(cfl)+'\n')
    f.write(str(testcase)+'\n')
    f.close()

    # run solver
    os.system(exe+" > out")
    os.system("tail -n1 out")
    os.system("tail -n1 out >> err.dat")

# read error from file
data = np.loadtxt('err.dat')
plt.autoscale(tight=False)
plt.loglog(data[:,2],data[:,4],'o-', linewidth=2)
plt.loglog(data[:,2],data[:,5],'*--',linewidth=2)
plt.loglog(data[:,2],data[:,2]**3,'-',linewidth=2)
plt.xlabel('h')
plt.ylabel('Error')
plt.legend(('L1','L2','Linf','slope_p94'))
plt.savefig('error_m2.pdf')
plt.show()
