%------------------- main -----------------------------------
% mgvee.m
function main_mgvee(problem)

% % V-cycle scheme to solve
% %     - Delta u = f on [0,1]
% % with boundary values given by the function g(x).
% N is the number of subintervals
% There are N+1 points
global N; N = 128; h = 1/N; tol = 1e-6;

% lmax determines the coarsest grid level; original grid = level 1
% e.g., lmax = 4 means one has to restrict
% to coarser grids 3 times
global lmax; lmax = 6;

% initial guess
j = 0:N; x = j*h;

if problem == 1
   initguess = sin(20*pi*x);
   vexact = sin(2*pi*x);    % exact solution
   f = 4*pi*pi*sin(2*pi*x); % right-hand side
elseif problem == 2
   initguess = 0*x;
   vexact = 0.5*x.*(1-x);
   f = ones(1,N+1);
end

v = initguess;

% main engine
ctr = 0; ctrmax = 50;
rfin_norm_init = norm(compresidual(v,f,1),2)/sqrt(length(f));
rfin_norm_old=rfin_norm_init; rfin_norm=rfin_norm_init;
f_norm = norm(f,2)/sqrt(length(f));
fprintf('Initial residual norm = %6.10d\n',rfin_norm_init)

while rfin_norm > tol * f_norm && ctr < ctrmax
   [vnew,rfin_norm] = vcycle (v, f, 1);
   conv_fact_r=rfin_norm/rfin_norm_old;
   rfin_norm_old=rfin_norm;
   fprintf('convergence factor based on residual is %6.10d\n',conv_fact_r);
   v = vnew;
   ctr = ctr + 1;
end
close all
plot(x,vexact,x,v); xlabel('x'); ylabel('v'); legend('Exact','FD'); grid on;
title(strcat('N = ',num2str(N)))
fprintf('The norm of the solution error is %6.10d\n', norm(v-vexact,2)/sqrt(length(v)))
fprintf('The norm of the solution residual is %6.10d\n', rfin_norm)
fprintf('The number of iterations required to satisfy tolerance is %d\n', ctr)