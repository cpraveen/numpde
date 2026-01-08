import numpy as np

def errornorm(h, u1, u2):
    error = u1 - u2
    return np.sqrt(h * np.sum(error[0:-1]**2))
