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

# Domain
xmin, xmax = 0.0, 2.0*pi

# Grid of n points
n = 20
h = (xmax - xmin)/(n - 1)
x = linspace(xmin,xmax,n)

# array for solution
u = zeros(n)

# BC for first and last points
u[0]  = uexact(x[0])
u[-1] = uexact(x[-1])

b = h**2 * f(x[1:-1])
b[0]  += u[0]
b[-1] += u[-1]
u[1:-1] = tdma(2*ones(n-2),-ones(n-2),-ones(n-2),b)
print("Max error = ", abs(uexact(x)-u).max())

# Exact solution on fine mesh for plotting
xe = linspace(xmin, xmax, 100); ue = uexact(xe)

# Plot exact and numerical solution
plt.plot(xe,ue,x,u,'o')
plt.legend(('Exact solution','Numerical solution'))
plt.xlabel('x')
plt.ylabel('u')
plt.show()
