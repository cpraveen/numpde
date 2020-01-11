"""
Solve u_t + u_x = 0 with periodic bc
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

# sawtooth
def sw1(x):
    return 0.5*(np.pi - x)

# truncated fourier series of sawtooth
def sw2(x):
    y = np.empty_like(x)
    for w in range(1,21):
        y += np.sin(w*x) / w
    return y

# truncated and smoothed fourier series of sawtooth
def sw3(x):
    y = np.empty_like(x)
    for w in range(1, 31):
        y += (1 - (w/30)**2) * np.sin(w*x) / w
    return y

# Following convention for periodicity
# Grid points = N+1
# h = (xmax - xmin)/(N+1)
# x[0] = xmin, x[N] = xmax-h
def solve(N, cfl, Tf, uinit):
    xmin, xmax = 0.0, 2.0*np.pi
    a = 1.0

    h = (xmax - xmin)/(N+1)
    dt = cfl * h / np.abs(a)
    print('dt, h = ',dt,h)

    x = np.linspace(xmin, xmax-h, N+1) # Note: last point does not reach xmax
    u2 = uinit(x)
    u4 = uinit(x)
    u6 = uinit(x)

    # First time step: forward euler
    v2 = u2 - a*dt*Q2(h,u2)
    v4 = u4 - a*dt*Q4(h,u4)
    v6 = u6 - a*dt*Q6(h,u6)
    t  = dt

    while t < Tf:
        # leap-frog
        w2 = u2 - 2.0*a*dt*Q2(h,v2)
        w4 = u4 - 2.0*a*dt*Q4(h,v4)
        w6 = u6 - 2.0*a*dt*Q6(h,v6)

        u2[:] = v2; v2[:] = w2
        u4[:] = v4; v4[:] = w4
        u6[:] = v6; v6[:] = w6
        t = t + dt

    plt.figure(figsize=(15,5))
    plt.subplot(131)
    plt.plot(x,w2); plt.title('Q2')
    plt.subplot(132)
    plt.plot(x,w4); plt.title('Q4')
    plt.subplot(133)
    plt.plot(x,w6); plt.title('Q6')
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of points', default=128)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.01)
parser.add_argument('-Tf', type=float, help='Final time', default=np.pi)
parser.add_argument('-ic', choices=('sw1', 'sw2', 'sw3'),
                    help='Init cond', default='sw1')
args = parser.parse_args()

# Run the solver
if args.ic == "sw1":
    solve(args.N, args.cfl, args.Tf, sw1)
elif args.ic == "sw2":
    solve(args.N, args.cfl, args.Tf, sw2)
elif args.ic == "sw3":
    solve(args.N, args.cfl, args.Tf, sw3)
