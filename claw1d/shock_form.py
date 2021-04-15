'''
Burgers shock formation
'''
import numpy as np
import matplotlib.pyplot as plt

x1, x2 = -0.5, 0.5
u1, u2 = 2.0, 1.0

tc = (x2 - x1)/(u1 - u2)
xc = x1 + u1*tc
s  = 0.5*(u1 + u2)

# Piecewise continuous
def uinit(x):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < x1:
            u[i] = u1
        elif y > x2:
            u[i] = u2
        else:
            u[i] = ((x2-y)*u1 + (y-x1)*u2)/(x2-x1)
    return u

def uexact1(x,t):
    y1 = x1 + u1*t
    y2 = x2 + u2*t
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < y1:
            u[i] = u1
        elif y > y2:
            u[i] = u2
        else:
            u[i] = ((y2-y)*u1 + (y-y1)*u2)/(y2-y1)
    return u

def uexact2(x,t):
    u = np.empty_like(x)
    for i,y in enumerate(x):
        if y < xc + s*(t - tc):
            u[i] = u1
        else:
            u[i] = u2
    return u

def uexact(x,t):
    if t < tc:
        return uexact1(x,t)
    else:
        return uexact2(x,t)

xmin, xmax = -1.0, 4.0
Tf = 2.0
nx = 500
nt = 101

x = np.linspace(xmin, xmax, nx)
u = uinit(x)

fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot(x, uinit(x), lw=2)
ax.set_xlabel('x'); ax.set_ylabel('u')
plt.title('t = 0.0'); plt.grid(True); 
plt.draw(); plt.pause(1.0)

ts = np.linspace(0.0, Tf, nt)

for t in ts:
    line.set_ydata(uexact(x,t))
    plt.title('t = '+str(round(t,3)))
    plt.draw(); plt.pause(0.1)

plt.show()

