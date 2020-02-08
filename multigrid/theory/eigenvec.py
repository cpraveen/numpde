from param import *
import numpy as np
import matplotlib.pyplot as plt

def get_eigen_vec(n):
    v = np.zeros((n+2,n+2))
    x = np.linspace(0,1,n+2)

    for k in range(1,n+1):
        v[k,:] = np.sin(k*np.pi*x)
    return x, v

def plot_eigen_vec():
    n = 11
    x, v = get_eigen_vec(n)

    fig = plt.figure(figsize=(8,8))

    fig.add_subplot(5,2,1) # k = 1
    plt.plot(x,v[1,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,2) # k = 2
    plt.plot(x,v[2,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,3) # k = 3
    plt.plot(x,v[3,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,4) # k = 4
    plt.plot(x,v[4,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,5) # k = 5
    plt.plot(x,v[5,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,6) # k = 7
    plt.plot(x,v[7,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,7) # k = 8
    plt.plot(x,v[8,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,8) # k = 9
    plt.plot(x,v[9,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,9) # k = 8
    plt.plot(x,v[10,:],'o-');
    fig.add_subplot(5,2,10) # k = 9
    plt.plot(x,v[11,:],'o-');

    plt.show()

def plot_figs(x,v):
    fig = plt.figure(figsize=(8,8))

    fig.add_subplot(5,2,1) # k = 1
    plt.plot(x,v[1,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,2) # k = 2
    plt.plot(x,v[2,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,3) # k = 3
    plt.plot(x,v[3,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,4) # k = 4
    plt.plot(x,v[4,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,5) # k = 5
    plt.plot(x,v[5,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,6) # k = 7
    plt.plot(x,v[7,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,7) # k = 8
    plt.plot(x,v[8,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,8) # k = 9
    plt.plot(x,v[9,:],'o-'); plt.xticks([])
    fig.add_subplot(5,2,9) # k = 8
    plt.plot(x,v[10,:],'o-');
    fig.add_subplot(5,2,10) # k = 9
    plt.plot(x,v[11,:],'o-');
    plt.show()
