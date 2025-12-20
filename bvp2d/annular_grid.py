import numpy as np
import matplotlib.pyplot as plt

m = 30 # radial points
n = 20 # angular points

R1, R2 = 1.0, 3.0

# Radial grid
r = np.linspace(R1,R2,m)
# Angular grid
theta = np.linspace(0.0,2.0*np.pi,n)
# Cartesian coordinates
x = np.outer(r, np.cos(theta))
y = np.outer(r, np.sin(theta))

# Plot angular lines
for i in range(m):
    plt.plot(x[i,:],y[i,:],'k-')
# Plot radial lines
for j in range(n):
    plt.plot(x[:,j],y[:,j],'k-')
plt.axis('equal')
plt.show()
