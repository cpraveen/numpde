import numpy as np

def const(x,t):
    if t < 1.0e-14:
        return np.ones(len(x))
    else:
        w = np.arange(1,100)
        w1 = 2*w - 1
        f  = np.kron(x,w1).reshape(len(x),len(w))
        u = np.sin(np.pi*f) * (np.exp(-w1**2 * np.pi**2 * t)/w1)
        return (4/np.pi) * np.sum(u,1)

# sin(pi*x)
def sine1(x,t):
    return np.exp(-np.pi**2*t) * np.sin(np.pi*x)

# sin(pi*x) + sin(4*pi*x)
def sine2(x,t):
    return (np.exp(-np.pi**2*t) * np.sin(np.pi*x)
            + np.exp(-16*np.pi**2*t) * np.sin(4*np.pi*x))

# Triangular hat function
def tri(x,t):
    u = np.empty_like(x)
    for i in range(len(x)):
        if x[i] < 0.5:
            u[i] = x[i]
        else:
            u[i] = 1.0 - x[i]
    return u

# Square hat function
def square(x,t):
    u = np.empty_like(x)
    for i in range(len(x)):
        if 0.25 < x[i] < 0.75:
            u[i] = 1.0
        else:
            u[i] = 0.0
    return u
