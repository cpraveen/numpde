%------------------------------------------------------------------------------
% Solve linear convection equation with periodic BC flux limiter scheme
% function lin_hyp_1d_periodic_fluxlimiter(N, cfl, lim, prob)
% N      = number of grid points
% cfl    = cfl number
% lim    = 0 first order upwind
%        = 1 Lax-wendroff
%        = 2 minmod
%        = 3 vanleer
%        = 4 superbee
% prob   = 1 sine
%        = 2 hat
%------------------------------------------------------------------------------
function lin_hyp_1d_periodic_fluxlimiter(N, cfl, lim, prob)

global a nu lambda limiter problem h

limiter = lim;
problem = prob;

a    = 1;

% Set domain and final time
if problem==1
   xmin = 0; xmax = 1; Tf   = 5;
elseif problem==2
   xmin = 0; xmax = 1; Tf   = 1;
else
   fprintf(1, 'Unknown problem %d\n', problem);
   pause
end

h      = (xmax - xmin)/(N-1);
dt     = cfl * h / abs(a);
nu     = a * dt / h;
lambda = dt/h;

fprintf(1,'N      = %d\n', N);
fprintf(1,'cfl    = %f\n', cfl);
fprintf(1,'limiter= %d\n', limiter);
fprintf(1,'h      = %f\n', h);
fprintf(1,'dt     = %f\n', dt);
fprintf(1,'nu     = %f\n', nu);

% Make grid
x  = linspace(xmin, xmax, N);

% Set initial condition
u = exact_solution(0, x);

t = 0;
while t < Tf
   u = update(u);
   t = t + dt;
   fe = exact_solution(t, x - a * t);
   plot(x, u, 'r-', x, fe, 'b--', 'LineWidth', 2)
   legend('Numerical', 'Exact')
   pause(0.001);
end

%------------------------------------------------------------------------------
% Exact solution
%------------------------------------------------------------------------------
function f = exact_solution(t,x)

global problem

if problem==1
   f = sin(2*pi*x);
elseif problem==2
   xx = x - floor(x);
   f  =  heaviside(xx-0.25)-heaviside(xx-0.75);
else
   fprintf(1, 'Unknown problem - %d\n', problem);
   pause
end

%------------------------------------------------------------------------------
% Limiter function
%------------------------------------------------------------------------------
function phi = LIMITER(dul, dur)

global limiter h

if limiter==0
   % first order upwind
   phi=0;
   return;
elseif limiter==1
   phi=1;
   % lax-wendroff
   return;
end

% Second order limited schemes
if dul*dur <= 0
   phi = 0.0;
else
   r = dul/dur;
   if limiter==2
      % minmod
      phi = min(r,1);
   elseif limiter==3
      % vanleer
      phi = 2*r/(1 + r);
   elseif limiter==4
      % superbee
      phi = max( min(2*r,1), min(r,2) );
   else
      fprintf(1, 'Unknown limiter %d\n', limiter);
      pause
   end
end

%------------------------------------------------------------------------------
% Numerical flux
%------------------------------------------------------------------------------
function f = num_flux(ujm1, uj, ujp1)

global a nu

dul = uj - ujm1;
dur = ujp1 - uj;

phi = LIMITER(dul, dur);

f = a * uj + 0.5 * a * (1 - nu) * phi * dur;

%------------------------------------------------------------------------------
% Backward difference in space
%------------------------------------------------------------------------------
function u = update(u)

global lambda

uold = u;
N = length(u);

res = zeros(size(u));

% flux across first face
uj = u(1); ujp1 = u(2); ujm1 = u(N-1);
flux = num_flux(ujm1, uj, ujp1);
res(2) = res(2) - flux;
res(N) = res(N) + flux;

for j=2:N-1
   uj = u(j); ujp1 = u(j+1); ujm1 = u(j-1);
   flux = num_flux(ujm1, uj, ujp1);
   res(j) = res(j) + flux;
   res(j+1) = res(j+1) - flux;
end

u(2:N) = u(2:N) - lambda * res(2:N);
u(1)   = u(N);
