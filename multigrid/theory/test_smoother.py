from param import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
from smoother import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-scheme', choices=('jacobi','gs'), help='Scheme',
                    default='jacobi')
args = parser.parse_args()

# grid of n+2 points
n = 64
h = 1.0/(n+1)
x = np.linspace(0.0,1.0, n+2) # grid
f = np.zeros(n+2) # rhs
# initial condition
v = (1.0/3.0)*(np.sin(np.pi*x) + np.sin(16*np.pi*x) + np.sin(32*np.pi*x))
fig = plt.figure()
ax = fig.add_subplot(111)
line1, = ax.plot(x, v, 'r')
line2, = ax.plot(x, v, 'b')
ax.set_xlabel('x')
ax.set_ylabel('v')
plt.grid(True)
plt.legend(('Initial', 'Current'))
plt.title('Scheme='+args.scheme+', n='+str(n))
plt.draw(); plt.pause(0.1)
wait = input("Press enter to continue ")

# for weighted jacobi
omega = 2.0/3.0

for i in range(10):
    if args.scheme == 'jacobi':
        v = wjacobi(h,v,f,omega,1)
    elif args.scheme == 'gs':
        v = gs(h,v,f,1)
    else:
        print('Unknown scheme')
        exit()
    line2.set_ydata(v)
    plt.draw(); plt.pause(1)

plt.show()
