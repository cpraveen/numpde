%------------------- main -----------------------------------
% mgvee.m

clear

% % V-cycle scheme to solve
% %     - Delta u = f on [0,1]
% % with boundary values given by the function g(x).

% N is the number of subintervals

global N; N = 128; h = 1/N; tol = 10e-6;

% lmax determines the coarsest grid level; original grid = level 1
% e.g., lmax = 4 means one has to restrict
% to coarser grids 3 times
global lmax; lmax = 6;


% initial guess
j = 0:N; initguess = sin(20*pi*j*h); v = initguess;

%zeros(size(initguess));

% exact solution
vexact = sin(2*pi*j*h);

% right-hand side
f = 4*pi*pi*sin(2*pi*j*h);


% main engine
relaterr = 10; ctr = 0; rfin_norm_old=1.0; rfin_norm=1.0;

relaterr_old=1.0;


while rfin_norm > tol
    [vnew,rfin_norm] = vcycle (v, f, 1);


% some output info
    relaterr = norm(vexact-vnew,2)/(norm(vexact,2)+.1);
    fprintf('relative error is %6.10d\n',relaterr);

    conv_fact_r=rfin_norm/rfin_norm_old;
    rfin_norm_old=rfin_norm;
    fprintf('convergence factor based on residual is %6.10d\n',conv_fact_r);

    v = vnew;
    ctr = ctr +1;
end

fprintf('The norm of the solution error is %6.10d\n', norm(v-vexact,2))
fprintf('The norm of the solution residual is %6.10d\n', rfin_norm)
fprintf('The number of iterations required to satisfy tolerance is %d\n', ctr)

