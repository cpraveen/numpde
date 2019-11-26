import numpy as np

def smooth(x):
    return np.sin(2*np.pi*x)

def hat(x):
    u = np.empty_like(x)
    for i,xx in enumerate(x):
        xx = np.abs(xx - int(xx)) # xx in [0,1]
        if np.abs(xx-0.5) < 0.25:
            u[i] = 1.0
        else:
            u[i] = 0.0
    return u

def sine(x):
    return 1.0 + np.sin(2*np.pi*x) + np.sin(20*np.pi*x)

