import numpy as np
import matplotlib.pyplot as plt


N = 500

zr = np.linspace(-2.5,1.0,N)
zi = np.linspace(-2.0,2.0,N)

zr,zi = np.meshgrid(zr,zi)
z = zr + 1j * zi

G = 1.0 + z + z**2/2.0

plt.figure(figsize=(6,5))
cs = plt.contourf(zr,zi,np.abs(G),levels=30,cmap='jet')
plt.colorbar(cs)
plt.contour(zr,zi,np.abs(G),levels=[0.9999,1,1.0001],linewidths=2)
plt.plot([0,0],[-1.73,1.73],'sg')
plt.axhline(0, color='white')
plt.axvline(0, color='white')
plt.text(-1.5,0,'Stable region',verticalalignment='center',
         color='white',bbox=dict(facecolor='grey'))
plt.text(-4, 0, 'Unstable region', verticalalignment='center',
         color='white', rotation=90,bbox=dict(facecolor='red'))
plt.title('RK3: Magnitude of amplification factor')
plt.xlabel('Real z')
plt.ylabel('Imaginary z')
plt.show()
