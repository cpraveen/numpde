"""
Convergence of finite difference for 
u(x) = sin(2*pi*x), x in [0,1]
First and second derivative using central difference
"""
import numpy as np
import matplotlib.pyplot as plt

xmin, xmax = 0.0, 1.0

def fun(x):
    return np.sin(2.0*np.pi*x)
    
def dfun(x):
    return 2.0*np.pi*np.cos(2.0*np.pi*x)

def ddfun(x):
    return -(2.0*np.pi)**2 * np.sin(2.0*np.pi*x)

npt = [50,100,200,400,800]
err1 = []; err2 = []; hh=[]
for n in npt:
    x = np.linspace(xmin,xmax,n)
    h = (xmax - xmin)/(n-1); hh.append(h)
    f = fun(x)
    # first derivative
    df = (f[2:] - f[0:-2])/(2*h)
    dfe= dfun(x)
    err1.append(np.abs(dfe[1:-1] - df).max())
    # second derivative
    ddf = (f[2:] - 2*f[1:-1] + f[0:-2])/(h*h)
    ddfe= ddfun(x)
    err2.append(np.abs(ddfe[1:-1] - ddf).max())
    
plt.loglog(hh,err1,'o-',hh,err2,'*-')
plt.xlabel('h')
plt.ylabel('Maximum error')
plt.legend(('First derivative','Second derivative'),loc='lower right')
plt.show()
