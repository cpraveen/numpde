"""
Solve u_t + f(u)_x = 0  for f(u) = u^2/2
Schemes are written in finite difference form
"""
import numpy as np
import matplotlib.pyplot as plt
import argparse

ul, ur = 1.2, 0.4
s      = 0.5*(ul + ur)

def uexact(t,x):
    u = np.zeros(len(x))
    for i in range(len(x)):
        if x[i] < 0.25 + s*t:
            u[i] = ul
        else:
            u[i] = ur
    return u

def update_fdup(lam, u):
    unew = np.empty_like(u)
    unew[0] = u[0];
    for i in range(1,len(u)-1):
        if u[i] >= 0.0:
            fx = u[i] * (u[i] - u[i-1])
        else:
            fx = u[i] * (u[i+1] - u[i])
        unew[i] = u[i] - lam*fx
    unew[-1] = u[-1];
    return unew

def update_lw(lam, u):
    unew = np.empty_like(u)
    f = 0.5*u*u
    unew[0] = u[0];
    for i in range(1,len(u)-1):
        al, ar = 0.5*(u[i-1]+u[i]), 0.5*(u[i+1]+u[i])
        unew[i] = u[i] - 0.5*lam*(f[i+1]-f[i-1]) \
                  + 0.5*lam**2*(ar*(u[i+1]-u[i]) - al*(u[i]-u[i-1]))
    unew[-1] = u[-1];
    return unew

def update_lf(lam, u):
    f = 0.5*u*u
    unew = np.empty_like(u)
    unew[0] = u[0];
    for i in range(1,len(u)-1):
        unew[i] = 0.5*(u[i-1]+u[i+1]) - 0.5*lam*(f[i+1] - f[i-1])
    unew[-1] = u[-1];
    return unew

def solve(N, cfl, scheme, Tf):
    xmin, xmax = 0.0, 1.0

    x = np.linspace(xmin, xmax, N)
    h = (xmax - xmin)/(N-1)
    u = uexact(0.0, x)
    dt= cfl * h / np.max(u)
    lam = dt/h

    fig = plt.figure()
    ax = fig.add_subplot(111)
    line1, = ax.plot(x, u, 'o')
    line2, = ax.plot(x, u, 'r')
    ax.set_xlabel('x'); ax.set_ylabel('u')
    plt.legend(('Numerical','Exact'))
    plt.title('N='+str(N)+', CFL='+str(cfl)+', Scheme='+scheme)
    plt.axis([0.0, 1.0, 0.0, 1.4])
    plt.draw(); plt.pause(0.1)
    wait = input("Press enter to continue ")

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='FDUP':
            u = update_fdup(lam, u)
        elif scheme=='LF':
            u = update_lf(lam, u)
        elif scheme=='LW':
            u = update_lw(lam, u)
        else:
            print("Unknown scheme: ", scheme)
            return
        t += dt; it += 1
        line1.set_ydata(u)
        line2.set_ydata(uexact(t,x))
        plt.draw(); plt.pause(0.1)
    plt.show()

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', type=int, help='Number of cells', default=100)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('FDUP','LF','LW'), help='Scheme', default='FDUP')
parser.add_argument('-Tf', type=float, help='Final time', default=0.5)
args = parser.parse_args()

# Run the solver
solve(args.N, args.cfl, args.scheme, args.Tf)
