"""
Convergence of finite difference for first derivative
"""
import numpy as np
import matplotlib.pyplot as plt

xmin, xmax = -1.0, 1.0

def fun(x):
    return x + np.exp(-np.abs(x)**3)
    
def dfun(x):
    return 1.0 - 3.0 * np.exp(-np.abs(x)**3) * x**2 * np.sign(x)

npt = [50,100,200,400,800]
errb = []; errf = []; errc = []; hh=[]
for n in npt:
    x = np.linspace(xmin,xmax,n)
    h = (xmax - xmin)/(n-1); hh.append(h)
    f = fun(x)
    # backward difference
    dfb = (f[1:-1] - f[0:-2])/h
    # forward difference
    dff = (f[2:] - f[1:-1])/h
    # central difference
    dfc = (f[2:] - f[0:-2])/(2*h)
    # exact derivative
    dfe= dfun(x[1:-1])
    # errors
    #plt.semilogy(x[1:-1], np.abs(dfe-dfc),'-o')
    errb.append(np.abs(dfe - dfb).max())
    errf.append(np.abs(dfe - dff).max())
    errc.append(np.abs(dfe - dfc).max())
    
plt.loglog(hh,errb,'o-',hh,errf,'*-',hh,errc,'s--')
plt.xlabel('h')
plt.ylabel('Maximum error')
plt.legend(('Backward','Forward','Central'),loc='lower right')
plt.title('First derivative for $u(x)=|x|^3$, $x\in[-1,1]$')
plt.show()
