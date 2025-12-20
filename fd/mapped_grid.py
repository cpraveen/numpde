import numpy as np
import matplotlib.pyplot as plt

f = lambda x: 0.5*(1 - np.cos(np.pi*x))
df= lambda x: 0.5*np.pi*np.sin(np.pi*x)

plt.figure()
xi = np.linspace(0,1,100)
plt.plot(xi,f(xi),label='f')
plt.plot(xi,df(xi),label='df')
plt.xlabel('$\\xi$')
plt.legend()

plt.figure()
n = 20
xi = np.linspace(0.0,1.0,n)
x = f(xi)
plt.plot(xi,0*xi,'s',label='$\\xi$ space')
plt.plot(x,1+0*x,'o',label='$x$ space')
plt.legend()
plt.title('Grid with '+str(n)+' points')
plt.yticks([])
plt.show()
