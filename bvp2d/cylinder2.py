'''
Generates polar mesh with nearly isotropic cells
Similar to cylinder.py but uses a map from [0,1] --> [R1,R2]
'''
import numpy as np
import matplotlib.pyplot as plt

n = 100 # no. angular points

R1, R2 = 1.0, 50.0

m = (n/(2*np.pi)) * np.log(R2/R1)
m = int(m)

# Map uniform partition of xi in [0,1] to r in [R1,R2]
xi = np.linspace(0.0,1.0,m)
r = R1 * np.exp(xi * np.log(R2/R1))
print("outer radius = ", r[-1])

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
