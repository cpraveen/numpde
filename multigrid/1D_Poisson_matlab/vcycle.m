%------------------------------------------------------
% vcycle.m


function [vcycleout,rfin_norm] = vcycle (v, f, L)
    global N;
    global lmax;

    % number of iterations in one level
    numiter = 2;
    v = wjacobi(v, f, numiter, L);

    % if not in yet the coarsest grid, restrict
    % if already in coarsest grid, relax and leave
    if L ~= lmax
  % should be rh instead of r2h, it is still fine grid
        r2h = compresidual (v,f,L);
        f2h = restrictfw (r2h,L);   %output is now at level L+1
  % zero initial guess
        v2h = zeros(size(f2h));
  % recursion:
        v2h_new = vcycle (v2h, f2h, L+1);
    else
        v = wjacobi (v, f, numiter, L);
        vcycleout = v; return
    end
  % prolongate error v2h_new
        errh = prolongate(v2h_new-v2h, L+1);    %output is now at level L
        v = v + errh;  %% error correction

    % relax and then leave
        v = wjacobi (v, f, numiter, L);
        vcycleout = v;

        rfin_norm = norm(compresidual(v,f,L),2);

% end
