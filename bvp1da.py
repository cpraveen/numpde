"""
Solve -u'' = f using finite difference
"""
from numpy import pi, sin, linspace, zeros, ones, abs
import matplotlib.pyplot as plt

# RHS function
def f(x):
    return sin(x)

# exact solution
def uexact(x):
    return sin(x)

# Solve using Thomas Tridiagonal algorithm
# A x = f
def solve(f):
    n = len(f)
    # construct LU decomposition
    a = 2.0*ones(n)
    b =-1.0*ones(n)
    c =-1.0*ones(n)
    for i in range(1,n):
        b[i] = b[i]/a[i-1]
        a[i] = a[i] - b[i]*c[i-1]
    # solution begins
    x = zeros(n)
    # solve L y = f
    x[0] = f[0]
    for i in range(1,n):
        x[i] = f[i] - b[i]*x[i-1]
    # solve U x = y
    x[-1] = x[-1]/a[-1]
    for i in range(n-2,-1,-1):
        x[i] = (x[i] - c[i]*x[i+1])/a[i]
    return x

xmin, xmax = 0.0, 2.0*pi
n = 20
h = (xmax - xmin)/(n - 1)

x = linspace(xmin,xmax,n)
u = zeros(n)
# BC for first and last points
u[0]  = uexact(x[0])
u[-1] = uexact(x[-1])

b = h**2 * f(x[1:-1])
b[0]  += u[0]
b[-1] += u[-1]
u[1:-1] = solve(b)
print "Max error = ", abs(uexact(x)-u).max()

xe = linspace(xmin, xmax, 100); ue = uexact(xe)
plt.plot(xe,ue,x,u,'o-')
plt.legend(('Exact solution','Numerical solution'))
plt.xlabel('x')
plt.ylabel('u')
plt.show()
