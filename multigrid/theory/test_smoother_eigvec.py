from param import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
from eigenvec import *
from smoother import *

# Get arguments
parser = argparse.ArgumentParser()
parser.add_argument('-scheme', choices=('jacobi','gs'), help='Scheme',
                    default='jacobi')
args = parser.parse_args()

n = 11
x, v = get_eigen_vec(n)
plot_figs(x,v)

h = 1.0/(n+1)
omega = 2.0/3.0
f = np.zeros(n+2)

for i in range(10):
    for k in range(1,n+1):
        if args.scheme == 'jacobi':
            v[k,:] = wjacobi(h,v[k,:],f,omega,1)
        elif args.scheme == 'gs':
            v[k,:] = gs(h,v[k,:],f,1)
        else:
            print('Unknown scheme')
            exit()
    plot_figs(x,v)
