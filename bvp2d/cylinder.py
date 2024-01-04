'''
Generates polar mesh with nearly isotropic cells
'''
import numpy as np
import matplotlib.pyplot as plt

n = 100 # no. angular points

R1, R2 = 1.0, 50.0

# Radial grid
r = np.array([R1])

while r[-1] < R2:
    dr = 2.0 * np.pi * r[-1] / n
    r = np.append(r, r[-1] + dr)

# no. of radial points
m = len(r)

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
