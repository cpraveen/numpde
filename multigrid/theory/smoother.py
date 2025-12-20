import numpy as np

# weighted Jacobi itertions
# omega=1 is standard Jacobi
# omega=2/3 is optimal
def wjacobi(h, v, f, omega, niter):
    n = len(v) - 2
    for i in range(niter):
        vtemp = v.copy()
        v[1:n+1] = 0.5*(v[0:n] + v[2:]) + 0.5 * h**2 * f[1:n+1]
        v = omega * v + (1.0-omega) * vtemp
    return v

# Gauss-Seidel scheme
def gs(h, v, f, niter):
    n = len(v) - 2
    for i in range(niter):
        for j in range(1,n+1): # TODO: loop is slow
            v[j] = 0.5*(v[j-1] + v[j+1]) + 0.5 * h**2 * f[j]
    return v
