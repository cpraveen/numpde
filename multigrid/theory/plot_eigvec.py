from param import *
import numpy as np
import matplotlib.pyplot as plt
import argparse
from eigenvec import *


# grid of n+2 points
n = 11
x, v = get_eigen_vec(n)
plot_figs(x, v)
plt.show()