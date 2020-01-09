%------------------------------------------------------------------------------
% Solve linear convection equation with periodic BC
% function linhyp1d_periodic_central(N, cfl, rks)
% N      = number of grid points
% cfl    = cfl number
% rks    = 1, 2, 3, 4
%------------------------------------------------------------------------------
function linhyp1d_periodic_central(N, cfl, rks)

xmin = 0;
xmax = 1;
a    = 1;
Tf   = 10;

h  = (xmax - xmin)/(N-1);
dt = cfl * h / abs(a);
nu = a * dt / h;

fprintf(1,'N      = %d\n', N);
fprintf(1,'cfl    = %f\n', cfl);
fprintf(1,'RK stag= %d\n', rks);
fprintf(1,'h      = %f\n', h);
fprintf(1,'dt     = %f\n', dt);
fprintf(1,'nu     = %f\n', nu);

% Make grid
x  = linspace(xmin, xmax, N);

% Initial condition
f = @(x) sin(2*pi*x);

% Set initial condition
u = f(x);

t = 0;
while t < Tf
   uold = u;
   for rk=1:rks
      fact = nu/(rks-rk+1);
      res = residual(u);
      u = uold - fact*res;
   end
   t = t + dt;
   fe = f(x - a * t);
   ti=strcat('t = ',num2str(t));
   plot(x, u, 'r-', x, fe, 'b--', 'LineWidth', 2)
   grid on
   legend('Numerical', 'Exact')
   title(ti)
   pause(0.1);
end

%------------------------------------------------------------------------------
% Backward difference in space
%------------------------------------------------------------------------------
function res = residual(u)

N = length(u);
res = zeros(1,N);

res(1)     = 0.5*(u(2) - u(N-1));
res(2:N-1) = 0.5*(u(3:N) - u(1:N-2));
res(N)     = 0.5*(u(2) - u(N-1));
