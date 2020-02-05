% Solve 2-D poisson equation using direct solver
n=25;
xmin=0; xmax=1; ymin=xmin; ymax=xmax;
h = (xmax-xmin)/(n-1);
m=n-2;
I = speye(m); e = ones(m,1);
D = spdiags([e -2*e e],[-1 0 1],m,m)/h^2;
A = -(kron(I,D) + kron(D,I));
x=linspace(xmin,xmax,n);
y=linspace(ymin,ymax,n);
[X,Y]=meshgrid(x,y);
f=2*(2*pi)^2*sin(2*pi*X).*sin(2*pi*Y);
utmp = A \ reshape(f(2:end-1,2:end-1),[m*m,1]);
u = zeros(n,n);
u(2:end-1,2:end-1) = reshape(utmp,[m,m]);
figure(1); contourf(X,Y,u,25); title('Numerical solution'); colorbar;
% Exact solution
ue=sin(2*pi*X).*sin(2*pi*Y);
figure(2); contourf(X,Y,u-ue,25); title('Error'); colorbar;
