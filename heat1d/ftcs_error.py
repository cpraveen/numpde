import numpy as np
import matplotlib.pyplot as plt

h   = 1.0/50
lam = 0.45
p   = 1
xi  = np.pi * p * h
a   = 1.0 - 4.0*lam*np.sin(0.5*xi)**2
b   = lam * h**2 * np.pi**2
print(a,b)

error = lambda n: np.abs(a**n - np.exp(-b*n))

n = np.arange(0,2000)
plt.plot(n, error(n))
plt.xlabel('n')
plt.ylabel('Error')
plt.title('h = '+str(h)+', lam = '+str(lam)+', p = '+str(p))
plt.show()
