"""
Solve u_t + u_x = 0 with periodic bc
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

def f(x):
    return np.sin(2*np.pi*x)

def update_ftbs(nu, u):
    unew = np.empty_like(u)
    for i in range(1,len(u)):
        unew[i] = (1-nu)*u[i] + nu*u[i-1]
    unew[0] = unew[-1]
    return unew

def update_ftfs(nu, u):
    unew = np.empty_like(u)
    for i in range(0,len(u)-1):
        unew[i] = (1+nu)*u[i] - nu*u[i+1]
    unew[-1] = unew[0]
    return unew

def update_ftcs(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] + 0.5*nu*(u[-2] - u[1])
    for i in range(1,len(u)-1):
        unew[i] = u[i] + 0.5*nu*(u[i-1] - u[i+1])
    unew[-1] = unew[0]
    return unew

def update_lf(nu, u):
    unew = np.empty_like(u)
    unew[0] = 0.5*(u[-1] + u[1]) + 0.5*nu*(u[-2] - u[1])
    for i in range(1,len(u)-1):
        unew[i] = 0.5*(u[i-1] + u[i+1]) + 0.5*nu*(u[i-1] - u[i+1])
    unew[-1] = unew[0]
    return unew

def update_lw(nu, u):
    unew = np.empty_like(u)
    unew[0] = u[0] - 0.5*nu*(u[1]-u[-2]) + 0.5*nu**2*(u[-2]-2*u[0]+u[1])
    for i in range(1,len(u)-1):
        unew[i] = u[i] - 0.5*nu*(u[i+1]-u[i-1]) + 0.5*nu**2*(u[i-1]-2*u[i]+u[i+1])
    unew[-1] = unew[0]
    return unew

def solve(N, cfl, scheme, Tf):
    xmin, xmax = 0.0, 1.0
    a          = 1.0

    h = (xmax - xmin)/N
    dt= cfl * h / np.abs(a)
    nu= a * dt / h

    x = np.linspace(xmin, xmax, N+1)
    u = f(x)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'r')
    line2, = ax.plot(x, u, 'b')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.draw(); plt.pause(0.1)
    wait = raw_input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='FTBS':
            u = update_ftbs(nu, u)
        elif scheme=='FTFS':
            u = update_ftfs(nu, u)
        elif scheme=='FTCS':
            u = update_ftcs(nu, u)
        elif scheme=='LF':
            u = update_lf(nu, u)
        elif scheme=='LW':
            u = update_lw(nu, u)
        else:
            print "Unknown scheme: ", scheme
            return
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(f(x-a*t))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('FTBS','FTFS','FTCS','LF','LW'), help='Scheme', default='FTBS')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.scheme, args.Tf)
