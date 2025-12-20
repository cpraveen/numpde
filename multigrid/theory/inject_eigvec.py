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
h = 1/(n+1)
xh, vh = get_eigen_vec(n)
plt.subplot(211)
plt.plot(xh, vh[5,:],'o-')
plt.plot(xh, -1.1*np.ones(n+2), 's-')
plt.grid(True)
plt.xticks([])
plt.title('$n=11, w_5^h$')
plt.ylim([-1.5,1.2])
for i in range(n+2):
    plt.text(i*h-0.01,-1.45,str(i))

n = 5
H = 1/(n+1)
xH, vH = get_eigen_vec(n)
plt.subplot(212)
plt.plot(xH, vH[5,:],'o-')
plt.plot(xH, -1.1*np.ones(n+2), 's-')
plt.grid(True)
plt.xlabel('x')
plt.title('$n=5, w_5^{2h}$')
plt.ylim([-1.5,1.2])
for i in range(n+2):
    plt.text(i*H-0.01,-1.45,str(i))

plt.show()

