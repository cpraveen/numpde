"""
Solve -u'' = f with dirichlet bc using finite difference
4th order monotone scheme of Bramble and Hubbard
https://doi.org/10.1090/S0025-5718-1964-0165702-X 
"""
from numpy import pi, sin, linspace, zeros, ones, abs, log
import matplotlib.pyplot as plt
from scipy.sparse import diags, csc_matrix
from scipy.sparse.linalg import spsolve

# RHS function
def f(x):
    return sin(x)

# exact solution
def uexact(x):
    return sin(x)

# Domain
xmin, xmax = 0.0, 2.0*pi

def error(n):
    # Grid of n points
    h = (xmax - xmin)/(n - 1)
    x = linspace(xmin,xmax,n)

    # array for solution
    u = zeros(n)

    # BC for first and last points
    u[0]  = uexact(x[0])
    u[-1] = uexact(x[-1])

    b = h**2 * f(x[1:-1])
    b[0]   += u[0]
    b[1 ]  -= u[0]/12.0
    b[-2]  -= u[-1]/12.0
    b[-1]  += u[-1]

    v0 = (5.0/2.0)*ones(n-2)
    v1 = -(4.0/3.0)*ones(n-3)
    v2 = (1.0/12.0)*ones(n-4)
    A = diags([v2,v1,v0,v1,v2], [-2,-1,0,1,2])
    A = A.tolil() # lil_matrix allows to change elements
    A[0,0] = 2.0; A[0,1] = -1.0; A[0,2] = 0.0
    A[-1,-3] = 0.0; A[-1,-2] = -1.0; A[-1,-1] = 2.0
    A = csc_matrix(A) # spsolve requires this format
    u[1:-1] = spsolve(A,b)
    return h, abs(uexact(x)-u).max()

# Compute error norm for different meshes
h = []; err = []
for n in [20,40,80,160,320]:
    h1, err1 = error(n)
    h.append(h1); err.append(err1)

# Compute convergence rate in L2 norm
print("h = %e   err = %e" % (h[0], err[0]))
for i in range(1,len(h)):
    p = log(err[i-1]/err[i])/log(2)
    print("h = %e   err = %e  rate = %f" % (h[i], err[i], p))

# Plot error norm vs h
plt.loglog(h, err, 'o-')
plt.xlabel('h')
plt.ylabel('Maximum Error')
plt.show()
plt.show()
