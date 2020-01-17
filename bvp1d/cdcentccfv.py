"""
Solve convection-diffusion eqn using central, cell-centered fvm
  -(1/Pe)*u'' + u' = 0 in (0,1)
   u(0) = 0,  u(1) = 1
"""
from numpy import zeros,linspace,exp
import matplotlib.pyplot as plt
import argparse
from tdma import *

parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', required=True)
parser.add_argument('-Pe', type=float, help='Peclet number', required=True)
args = parser.parse_args()
N, Pe = args.N, args.Pe

a, b = 0.0, 1.0

# Grid with N cells
h = (b-a)/N

# Mesh Peclet number
P = h * Pe

# Three diagonals
A,B,C = zeros(N),zeros(N),zeros(N)

A[0], C[0] = 1.0 + 6.0/P, 1.0 - 2.0/P
for i in range(1,N-1):
    B[i],A[i],C[i] = -(1.0 + 2.0/P), 4.0/P, 1.0 - 2.0/P
B[-1], A[-1] = -(1.0+2.0/P), -1.0 + 6.0/P

# Right hand side
F = zeros(N)
F[0]  = (2.0 + 4.0/P)*a
F[-1] = (-2.0 + 4.0/P)*b

# Solve
u = tdma(A,B,C,F)

# Cell centers of FV mesh
x = linspace(a+0.5*h, b-0.5*h, N)

# Exact solution on fine mesh for plotting
xe= linspace(a, b, 100)
ue= a + (b-a)*(exp((xe-1)*Pe) - exp(-Pe))/(1 - exp(-Pe))

# Plot numerical and exact solution
plt.plot(x,u,'o-',xe,ue)
plt.xlabel('x')
plt.ylabel('u')
plt.legend(('Numerical','Exact'),loc='upper left')
T='N='+str(N)+', Pe='+str(Pe)+', P='+str(P)
plt.title(T)
plt.show()
