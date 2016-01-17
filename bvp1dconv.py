"""
Solve -u'' = f using finite difference
"""
from numpy import pi, sin, linspace, zeros, ones, abs
import matplotlib.pyplot as plt
from tdma import *

# RHS function
def f(x):
    return sin(x)

# exact solution
def uexact(x):
    return sin(x)

# Solve the problem and compute error norm
def error(n):
    xmin, xmax = 0.0, 2.0*pi
    h = (xmax - xmin)/(n - 1)

    x = linspace(xmin,xmax,n)
    u = zeros(n)
    # BC for first and last points
    u[0]  = uexact(x[0])
    u[-1] = uexact(x[-1])

    b = h**2 * f(x[1:-1])
    b[0]  += u[0]
    b[-1] += u[-1]
    u[1:-1] = tdma(2*ones(n-2),-ones(n-2),-ones(n-2),b)
    return h, abs(uexact(x)-u).max()


h = []; err = []
for n in [20,40,80,160,320]:
    h1, err1 = error(n)
    print "h = %e   err = %e" % (h1, err1)
    h.append(h1); err.append(err1)

plt.loglog(h, err, 'o-')
plt.xlabel('h')
plt.ylabel('Maximum Error')
plt.show()
