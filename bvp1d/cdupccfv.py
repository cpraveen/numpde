"""
Solve convection-diffusion eqn using upwind, cell-centered fvm
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
h = (b-a)/N
P = h * Pe

# Construct diagonals of matrix
A,B,C = zeros(N),zeros(N),zeros(N)
A[0], C[0] = 1.0+2.0/P, -1.0/P
for i in range(1,N-1):
    B[i],A[i],C[i] = -(1.0+1.0/P), 1.0+2.0/P, -1.0/P
B[-1], A[-1] = -(1.0+1.0/P), 1.0+3.0/P

# construct rhs
F = zeros(N)
F[0]  = (1.0+1.0/P)*a
F[-1] = (2.0/P)*b

# solve
u = tdma(A,B,C,F)

# plot the solution
x = linspace(a+0.5*h, b-0.5*h, N)
xe= linspace(a, b, 100)
ue= a + (b-a)*(exp((xe-1)*Pe) - exp(-Pe))/(1 - exp(-Pe))

plt.plot(x,u,'o-',xe,ue)
plt.xlabel('x')
plt.ylabel('u')
plt.legend(('Numerical','Exact'),loc='upper left')
T='N='+str(N)+', Pe='+str(Pe)+', P='+str(P)
plt.title(T)
plt.show()
