'''
Lagrangian solution for Burger's equation
'''
import argparse
import numpy as np
import matplotlib.pyplot as plt

# Exponential
u1 = lambda x: 1.0 + np.exp(-10*x**2)

# Exponential
u2 = lambda x: 2.0 - np.exp(-10*x**2)

# Piecewise continuous
def u3(x):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < -0.5:
            u[i] = 2.0
        elif y > 0.5:
            u[i] = 1.0
        else:
            u[i] = - y + 1.5
    return u

# Piecewise continuous
def u4(x):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < -0.5:
            u[i] = 1.0
        elif y > 0.5:
            u[i] = 2.0
        else:
            u[i] = y + 1.5
    return u

# Jump
def u5(x):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < 0.0:
            u[i] = 2.0
        else:
            u[i] = 1.0
    return u

# Jump
def u6(x):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < 0.0:
            u[i] = 1.0
        else:
            u[i] = 2.0
    return u

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-nx', type=int, help='Number of x points', default=200)
parser.add_argument('-nt', type=int, help='Number of t points', default=50)
parser.add_argument('-Tf', type=float, help='Final time', default=1.0)
parser.add_argument('-ic', choices=('1','2','3','4','5','6'), 
                    help='Initial condition', default='1')
parser.add_argument('-pde', choices=('adv','burg'), 
                    help='PDE: lin adv or burger', default='burg')
args = parser.parse_args()

xmin, xmax = -1.0, 1.0
Tf = args.Tf
nx = args.nx
nt = args.nt

if   args.ic == '1':
    uinit = u1
elif args.ic == '2':
    uinit = u2
elif args.ic == '3':
    uinit = u3
elif args.ic == '4':
    uinit = u4
elif args.ic == '5':
    uinit = u5
elif args.ic == '6':
    uinit = u6

x = np.linspace(xmin, xmax, nx)
u = uinit(x)

ts = np.linspace(0.0, Tf, nt)

plt.figure()
for t in ts:
    plt.clf()
    if args.pde == 'adv':
        plt.plot(x + t, u, lw=2)
    else:
        plt.plot(x + u*t, u, lw=2)
    plt.xlabel('x'); plt.ylabel('u')
    plt.title('t = '+str(round(t,3)))
    plt.grid(True)
    plt.draw()
    if t == 0.0:
        plt.pause(2.0)
    else:
        plt.pause(0.1)

plt.show()

