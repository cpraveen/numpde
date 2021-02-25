"""
Solve using Thomas Tridiagonal algorithm
  A x = f
where A = diag(b,a,c)
"""
from numpy import zeros

# We modify the a,b arrays
# a = main diagonal
# b = sub diagonal
# c = super diagonal
def tdma(a,b,c,f):
    n = len(f)
    # construct LU decomposition
    for i in range(1,n):
        b[i] = b[i]/a[i-1]
        a[i] = a[i] - b[i]*c[i-1]
    # solution begins
    x = zeros(n)
    # solve L y = f; we use x to store the solution
    x[0] = f[0]
    for i in range(1,n):
        x[i] = f[i] - b[i]*x[i-1]
    # solve U x = y
    x[-1] = x[-1]/a[-1]
    for i in range(n-2,-1,-1):
        x[i] = (x[i] - c[i]*x[i+1])/a[i]
    return x
