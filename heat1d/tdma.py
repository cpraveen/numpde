"""
Solve using Thomas Tridiagonal algorithm
  A x = f
where A = diag(b,a,c)
"""
from numpy import zeros

# construct LU decomposition
def tdma1(a,b,c):
    n = len(a)
    for i in range(1,n):
        b[i] = b[i]/a[i-1]
        a[i] = a[i] - b[i]*c[i-1]

# solve
def tdma2(a,b,c,f):
    n = len(f)
    x = zeros(n)
    # solve L y = f
    x[0] = f[0]
    for i in range(1,n):
        x[i] = f[i] - b[i]*x[i-1]
    # solve U x = y
    x[-1] = x[-1]/a[-1]
    for i in range(n-2,-1,-1):
        x[i] = (x[i] - c[i]*x[i+1])/a[i]
    return x
