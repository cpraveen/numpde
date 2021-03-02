"""
Solve convection-diffusion eqn using central, cell-centered fvm
  -(1/Pe)*u'' + u' = 0 in (0,1)
   u(0) = 0,  u(1) = 1
"""
from numpy import zeros,linspace,exp,empty_like,floor
import matplotlib.pyplot as plt
import argparse
from tdma import *

# xf: location of cell faces
def solve(Pe,xf):
    N = len(xf) - 1 # number of cells
    xc = 0.5*(xf[0:-1] + xf[1:]) # cell centers
    hc = xf[1:] - xf[0:-1] # cell lengths
    hf = empty_like(xf) # h_{j+1/2}
    hf[0] = hc[0]
    hf[1:-1] = 0.5*(hc[0:-1] + hc[1:])
    hf[-1] = hc[-1]

    # Mesh Peclet number
    P = hf * Pe

    # Three diagonals
    A,B,C = zeros(N),zeros(N),zeros(N)

    A[0], C[0] = 0.5 + 2.0/P[0] + 1.0/P[1], 0.5 - 1.0/P[1]
    for i in range(1,N-1):
        B[i],A[i],C[i] = -(0.5+1.0/P[i]), (1.0/P[i]+1.0/P[i+1]), 0.5-1.0/P[i+1]
    B[-1], A[-1] = -(0.5+1.0/P[-2]), -0.5 + 2.0/P[-1] + 1.0/P[-2]

    a, b = 0.0, 1.0 # Dirichlet condition

    # Right hand side
    F = zeros(N)
    F[0]  = (1.0 + 2.0/P[0])*a
    F[-1] = (-1.0 + 2.0/P[-1])*b

    # Solve
    u = tdma(A,B,C,F)

    return xc, u

#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', required=True)
parser.add_argument('-Pe', type=float, help='Peclet number', required=True)
parser.add_argument('-grid', choices=('uniform','adaptive'), 
                    help='Grid type', required=True)
args = parser.parse_args()
N, Pe = args.N, args.Pe

a, b = 0.0, 1.0 # domain

# Make grid: N = number of cells
if args.grid == 'uniform':
    xf = linspace(a,b,N+1)
else:
    delta = 4.0/Pe  # Width of boundary layer
    x0, x1, x2 = a, b-delta, b
    N0 = int(floor(2*N/3))     # no. cells outside bl
    N1 = N - N0                # no. cells inside  bl

    xf = zeros(N+1) # FV face locations
    xf[0:N0+1] = linspace(x0,x1,N0+1) # N0 cells in [x0,x1]
    xf[N0:] = linspace(x1,x2,N1+1)    # N1 cells in [x1,x2]

x, u = solve(Pe,xf)

# Exact solution to estimate error
uex = a + (b-a)*(exp((x-1)*Pe) - exp(-Pe))/(1 - exp(-Pe))
print('Max error = ', abs(uex-u).max())

# Exact solution on fine mesh for plotting
xe= linspace(a, b, 100)
ue= a + (b-a)*(exp((xe-1)*Pe) - exp(-Pe))/(1 - exp(-Pe))

# Plot numerical and exact solution
plt.plot(x,u,'o-',label='Numerical')
plt.plot(xe,ue,label='Exact')
plt.plot(xf,0*xf,'k|',lw=2)
plt.xlabel('x')
plt.ylabel('u')
plt.legend(loc='upper left')
plt.grid(True)
T='N='+str(N)+', Pe='+str(Pe)
plt.title(T)
plt.show()
