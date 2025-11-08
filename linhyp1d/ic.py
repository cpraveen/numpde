import numpy as np

# Domain = [0,1]
def smooth(x):
    return np.sin(2*np.pi*x)

# Domain = [0,1]
def hat(x):
    u = np.empty_like(x)
    for i,xx in enumerate(x):
        xx = np.abs(xx - int(xx)) # xx in [0,1]
        if np.abs(xx-0.5) < 0.25:
            u[i] = 1.0
        else:
            u[i] = 0.0
    return u

# Domain = [0,1]
def sine(x):
    return 1.0 + np.sin(2*np.pi*x) + np.sin(20*np.pi*x)

# Domain = [-1,+1]
def mult(x):
    f = np.empty_like(x)
    for i,xx in enumerate(x):
        if xx > 1.0:
            y = xx - 2.0*int((xx+1.0)/2.0)
        elif xx < -1.0:
            y = xx + 2.0*int((1.0-xx)/2.0)
        else:
            y = xx
        if y > -0.8 and y < -0.6:
            f[i] = np.exp(-np.log(2.0)*(y+0.7)**2/0.0009);
        elif y > -0.4 and y < -0.2:
            f[i] = 1.0;
        elif y > 0.0 and y < 0.2:
            f[i] = 1.0 - np.abs(10.0*(y-0.1));
        elif y > 0.4 and y < 0.6:
            f[i] = np.sqrt(1.0 - 100.0*(y-0.5)**2);
        else:
            f[i] = 0.0;
    return f

# See Fig. 6.1, 6.2 of Leveque FVM book
# Domain [0,1]
def exphat(x):
    beta, x0 = 200.0, 0.3
    f = np.exp(-beta * (x-x0)**2) + (np.abs(x - 0.7) < 0.1) * 1.0
    return f

# Wave packet
# Domain: [-1,1]
def wpack(x):
    y = np.abs(x - np.floor(x)) # xx in [0,1]
    y = y - 0.5
    return np.sin(20 * np.pi * y) * np.exp(-50 * y**2)

# zero + noise
# Domain: [0,1]
def noise(x):
    y = 2.0 * np.random.rand(len(x)) - 1.0
    return 1.0e-16 * y
