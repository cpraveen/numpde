from param import *
import numpy as np
import matplotlib.pyplot as plt

n = 64
h = 1.0/(n+1)

plt.figure()
k = np.arange(1,n+1)
for omega in [1.0/3.0,2.0/3.0, 1.0]:
    lam = 1 - 2*omega*np.sin(0.5*np.pi*k*h)**2
    plt.plot(k/n,lam)

plt.plot([1/n,1],[1.0/3.0,1.0/3.0],'k--',lw=1)
plt.plot([1/n,1],[-1.0/3.0,-1.0/3.0],'k--',lw=1)
plt.legend(('$\omega=1/3$','$\omega=2/3$','$\omega=1$'))
plt.grid(True)
plt.xlabel('$k/n$')
plt.ylabel('$\lambda_k$')
plt.axis('tight')
plt.show()
