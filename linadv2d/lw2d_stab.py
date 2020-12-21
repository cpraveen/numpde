import numpy as np
from numpy import exp
import matplotlib.pyplot as plt

def amp_f(sigma_x,sigma_y,H1,H2):
    result = 1.0 - 0.5*sigma_x*(exp(1j*h1) - exp(-1j*h1)) \
                 - 0.5*sigma_y*(exp(1j*h2) - exp(-1j*h2)) \
                 + 0.5*sigma_x**2*(exp(-1j*h1) - 2.0 + exp(1j*h1)) \
                 + 0.5*sigma_y**2*(exp(-1j*h2) - 2.0 + exp(1j*h2)) \
                 + 0.25*sigma_x*sigma_y*(  exp(1j*h1)*exp(1j*h2) \
                                         - exp(1j*h1)*exp(-1j*h2) \
                                         - exp(-1j*h1)*exp(1j*h2) \
                                         + exp(-1j*h1)*exp(-1j*h2))
    return np.max(np.abs(result))

sigma_x_range = np.linspace(-1.0,1.0,200)
sigma_y_range = np.linspace(-1.0,1.0,200)
h1 = np.linspace(0,2*np.pi,100)
h2 = np.linspace(0,2*np.pi,100)
h1,h2 = np.meshgrid(h1,h2)
X, Y = [], []
A = [] # For all sigma_x for which (sigma_x,sigma_x) is a stable pair.
for sigma_x in sigma_x_range:
    for sigma_y in sigma_y_range:
        if amp_f(sigma_x,sigma_y,h1,h2) - 1.0 < 1e-5:
            X.append(sigma_x)
            Y.append(sigma_y)
            if np.abs(sigma_x-sigma_y) < 1e-4: # Checking
                A.append(sigma_x)
print('Highest sigma for which (sigma,sigma) is stable pair is ', np.max(A))
plt.scatter(X,Y,c='y',label='Stable Region')
plt.plot(A,A,label="$\sigma_x = \sigma_y$")
R = 2*0.36
Z = np.linspace(0.0,R,100) # Trying to fit the best line.
plt.plot(Z,R-Z,label='$\sigma_x+\sigma_y = 0.72$')
plt.title('Lax-Wendroff: Region of stable CFLs')
plt.xlabel('$\sigma_x$')
plt.ylabel('$\sigma_y$')
plt.axis('equal')
plt.legend()
plt.grid(True)
plt.show()
