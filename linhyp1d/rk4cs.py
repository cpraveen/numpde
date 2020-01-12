"""
Solve u_t + a u_x = 0 with periodic bc
central difference in space + RK4 in time
space order = 2, 4, 6
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse
from ic import *

# Backward difference
def Dm(h,u):
    return (u - np.roll(u,1))/h

# Forward difference
def Dp(h, u):
    return (np.roll(u,-1) - u)/h

# Central difference
# (u(i+1) - u(i-1))/(2 h)
def Q2(h,u):
    return (np.roll(u,-1) - np.roll(u,1))/(2.0*h)

# fourth order central difference
def Q4(h,u):
    a = Dm(h,u)
    b = Dp(h,a) # = (Dp*Dm)*u
    c = u - (h**2/6)*b
    return Q2(h,c)

# sixth order central difference
def Q6(h, u):
    a = Dm(h,u)
    b = Dp(h,a) # = (Dp*Dm)*u
    c = Dm(h,b)
    d = Dp(h,c) # = (Dp*Dm)^2 * u
    e = u - (h**2/6)*b + (h**4/30)*d
    return Q2(h,e)

# ic = sin(10 x)
def sin10(x):
    return np.sin(10*x)

# Following convention for periodicity
# Grid points = N+1
# h = (xmax - xmin)/(N+1)
# x[0] = xmin, x[N] = xmax-h
def solve(N, cfl, Tf, diff, uinit):
    xmin, xmax = 0.0, 2.0*np.pi
    a = 1.0

    h = (xmax - xmin)/(N+1)
    dt = cfl * h / np.abs(a)
    print('dt, h = ',dt,h)

    x = np.linspace(xmin, xmax-h, N+1) # Note: last point does not reach xmax
    u = uinit(x)

    err = np.array([[0,0]]) # collect (time, max error) in this array
    t = 0.0
    while t < Tf:
        k1 = -a * dt * diff(h,u)
        k2 = -a * dt * diff(h,u+0.5*k1)
        k3 = -a * dt * diff(h,u+0.5*k2)
        k4 = -a * dt * diff(h,u+k3)
        u += (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
        t = t + dt
        ue = uinit(x-a*t)
        err = np.append(err, [[t,np.abs(u-ue).max()]], axis=0)

    print('Final error =', err[-1,1])

    plt.figure()
    plt.plot(x,uinit(x-a*t),x,u)
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend(('Exact','Numerical'))
    plt.title('N='+str(N)+', cfl='+str(cfl))

    plt.figure()
    plt.plot(err[:,0],err[:,1])
    plt.xlabel('t')
    plt.ylabel('max error')
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of points', default=128)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.5)
parser.add_argument('-Tf', type=float, help='Final time', default=2*np.pi)
parser.add_argument('-scheme', type=int, choices=(2,4,6), help='derivative order', default=2)
args = parser.parse_args()

# select derivative scheme
if args.scheme == 2:
    diff = Q2
elif args.scheme == 4:
    diff = Q4
elif args.scheme == 6:
    diff = Q6
else:
    print('Unknown scheme=',args.scheme)

# Run the solver
solve(args.N, args.cfl, args.Tf, diff, sin10)
