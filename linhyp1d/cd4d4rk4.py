"""
Solve u_t + u_x = 0 with periodic bc and CD4 and 4th order dissipation
Compare to Fig. 7.1.3 in Gustafsson et al., Ed. 2
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

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

def rhs(a,h,eps,u):
    d = Dm(h,u)
    d = Dp(h,d)
    d = Dm(h,d)
    d = Dp(h,d)
    return -a*Q4(h,u) - eps * h**3 * d

def hat(x):
    return (np.abs(x - np.pi) < np.pi/3)*1.0 + 0.0

# Following convention for periodicity
# Grid points = N+1
# h = (xmax - xmin)/(N+1)
# x[0] = xmin, x[N] = xmax-h
def solve(N, cfl, eps, Tf, uinit):
    xmin, xmax = 0.0, 2.0*np.pi
    a = 1.0

    h = (xmax - xmin)/(N+1)
    dt = cfl * h / np.abs(a)
    print('dt, N, h = ',dt,N,h)

    x = np.linspace(xmin, xmax-h, N+1) # Note: last point does not reach xmax
    u = uinit(x)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_ylim(-0.5,1.5)
    line, = ax.plot(x, u, 'r')
    plt.grid(True)
    plt.title("N="+str(N+1))

    it,t = 0, 0.0
    while t < Tf:
        # RK4
        k0 = rhs(a, h, eps, u)
        k1 = rhs(a, h, eps, u+0.5*dt*k0)
        k2 = rhs(a, h, eps, u+0.5*dt*k1)
        k3 = rhs(a, h, eps, u+dt*k2)
        u += (dt/6)*(k0 + 2*k1 + 2*k2 + k3)
        t = t + dt
        it = it + 1
        if it%10 == 0:
            line.set_ydata(u)
            plt.pause(0.2)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of points', default=239)
parser.add_argument('-cfl', type=float, help='CFL number', default=2.0/3.0)
parser.add_argument('-eps', type=float, help='Epsilon', default=0.0)
parser.add_argument('-Tf', type=float, help='Final time', default=2*np.pi)
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.eps, args.Tf, hat)
