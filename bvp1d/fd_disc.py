"""
BVP with discontinuous coefficients using central FD
   - (a(x) u'(x))' = sin(2*pi*x)  in  (0,1)
   a(x) = al for x < 1/2
        = ar for x > 1/2
"""
from numpy import linspace,sin,pi,zeros,log
import matplotlib.pyplot as plt
from tdma import *

xmin,xmax = 0.0,1.0
al, ar = 1.0, 100.0

# viscosity coefficient
def afun(x):
    a = 0*x
    for i in range(len(x)):
        if x[i] <= 0.5:
            a[i] = al
        else:
            a[i] = ar
    return a

# right hand side function
def qfun(x):
    return sin(2.0*pi*x)

# right hand side function
def uexact(x):
    return sin(2.0*pi*x)/(4*pi**2*afun(x))

def solve(n):
    # we have n+1 grid points
    h = (xmax-xmin)/n
    x = linspace(xmin,xmax,n+1)
    a = afun(x)

    # Right hand side
    b = 2.0 * h**2 * qfun(x)
    # First and last equations are bc
    b[0]  = 0.0
    b[-1] = 0.0;

    A,B,C = zeros(n+1), zeros(n+1), zeros(n+1)

    A[0],C[0] = 1.0, 0.0
    for i in range(1,n):
        B[i] = -(a[i-1]+a[i])
        A[i] =  (a[i-1]+2*a[i]+a[i+1])
        C[i] = -(a[i]+a[i+1])
    B[-1],A[-1] = 0.0, 1.0

    u = tdma(A,B,C,b)

    xe = linspace(xmin,xmax,200)
    ue = uexact(xe)
    plt.plot(xe,ue,'-',x,u,'o--')
    plt.title('FD using '+str(n+1)+' points')
    plt.xlabel('x'); plt.ylabel('u')
    plt.legend(('Exact','FD'))
    print("Close plot window to continue")
    plt.show()
    return h, abs(u-uexact(x)).max()

# Compute error norm for different meshes
h = []; err = []
for n in [19,39,79,159,319]:
    h1, err1 = solve(n)
    h.append(h1); err.append(err1)

# Compute convergence rate in L2 norm
print("h = %e   err = %e" % (h[0], err[0]))
for i in range(1,len(h)):
    p = log(err[i-1]/err[i])/log(2)
    print("h = %e   err = %e  rate = %f" % (h[i], err[i], p))
