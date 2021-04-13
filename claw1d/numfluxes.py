"""
Solve u_t + f(u)_x = 0  for f(u) = u^2/2
Finite volume scheme
Smooth initial data
"""
import numpy as np

# Burgers flux
def flux(u):
    return 0.5*u*u

# Central flux
def flux_central(lam,u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        nf[i] = 0.5*(f[i-1] + f[i])
    return nf

# Lax-Friedrich flux
def flux_lf(lam,u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        nf[i] = 0.5*(f[i-1] + f[i]) - 0.5*(u[i] - u[i-1])/lam
    return nf

# Local Lax-Friedrich flux
def flux_llf(u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        a     = np.abs([u[i-1], u[i]]).max()
        nf[i] = 0.5*(f[i-1] + f[i]) - 0.5*a*(u[i] - u[i-1])
    return nf

# Lax-Wendroff flux
def flux_lw(lam,u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        a     = 0.5*(u[i-1]+u[i])
        nf[i] = 0.5*(f[i-1] + f[i]) - 0.5*lam*a*(f[i] - f[i-1])
    return nf

# Roe flux
def flux_roe(u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        a     = np.abs(0.5*(u[i-1]+u[i]))
        nf[i] = 0.5*(f[i-1] + f[i]) - 0.5*a*(u[i] - u[i-1])
    return nf

# Roe flux with entropy fic
def flux_eroe(u):
    n  = len(u)
    f  = flux(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        delta = 0.5*np.abs(u[i]-u[i-1]) if u[i-1] < u[i] else 0.0
        a     = np.abs(0.5*(u[i-1]+u[i]))
        a     = delta if a < delta else a
        nf[i] = 0.5*(f[i-1] + f[i]) - 0.5*a*(u[i] - u[i-1])
    return nf

# Godunov flux
def flux_god(u):
    n  = len(u)
    nf = np.zeros(n+1)
    for i in range(1,n):
        u1 = max(0.0, u[i-1])
        u2 = min(0.0, u[i])
        nf[i] = max(flux(u1), flux(u2))
    return nf
