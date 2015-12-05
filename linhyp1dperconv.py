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

    t, it = 0.0, 0
    while t < Tf:
        if scheme=='FTBS':
            u = update_ftbs(nu, u)
        elif scheme=='LW':
            u = update_lw(nu, u)
        else:
            print "Unknown scheme: ", scheme
            return
        t += dt; it += 1

    err = np.abs(u - f(x-a*t))
    em  = np.max(err)
    e1  = h*err[0] + h*np.sum(err[1:-2])
    err = err**2
    e2  = np.sqrt(h*err[0] + h*np.sum(err[1:-2]))
    return em,e1,e2

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-N', metavar='N', type=int, nargs='+', help='Number of cells', required=True)
parser.add_argument('-cfl', type=float, help='CFL number', default=0.9)
parser.add_argument('-scheme', choices=('FTBS','LW'), help='Scheme', default='FTBS')
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
args = parser.parse_args()

# Run the solver
emax = np.empty(len(args.N))
e1   = np.empty(len(args.N))
e2   = np.empty(len(args.N))
i    = 0
for N in args.N:
    print "Running for cells = ", N
    a,b,c = solve(N, args.cfl, args.scheme, args.Tf)
    emax[i] = a
    e1[i]   = b
    e2[i]   = c
    i += 1

plt.loglog(args.N, e1  , '*-')
plt.loglog(args.N, e2  , 's-')
plt.loglog(args.N, emax, 'o-')
plt.xlabel('N')
plt.ylabel('Error norm')
plt.title('Scheme='+args.scheme+', CFL='+str(args.cfl))
plt.legend(('$L_1$ norm','$L_2$ norm','$L_\infty$ norm'))
plt.show()
