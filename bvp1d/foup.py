"""
Solve using cc finite volume upwind scheme
   du/dx = q(x), x in (0,1)
   u(0) = a
"""
from numpy import linspace,zeros,sin,cos,pi
import matplotlib.pyplot as plt

xmin,xmax = 0.0, 1.0
a = 1.0
n = 20

# Right hand side function
def qfun(x):
    q = 0*x
    for i in range(len(x)):
        if x[i] < 0.5:
            q[i] = -1.0*sin(pi*x[i])
        else:
            q[i] = -2.0*sin(pi*x[i])
    return q

# Exact solution
def uexact(x):
    ue = 0*x
    for i in range(len(x)):
        if x[i] < 0.5:
            ue[i] = a + (cos(pi*x[i])-1)/pi
        else:
            ue[i] = a - 1.0/pi + 2.0*cos(pi*x[i])/pi
    return ue


h = (xmax-xmin)/n
x = linspace(xmin+0.5*h, xmax-0.5*h, n)
q = qfun(x)
u = zeros(n)

u[0] = a + h*q[0]
for i in range(1,n):
    u[i] = u[i-1] + h*q[i]

xe = linspace(xmin,xmax,100)
ue = uexact(xe)
plt.plot(x,u,'o',xe,ue,'-')
plt.xlabel('x')
plt.ylabel('u')
plt.legend(('Numerical','Exact'))
plt.title('Upwind scheme')
plt.show()
