%------------------------------------------------------
% compresidual.m

function residualout = compresidual (v, f, L)
   global N;
   n = N / 2^(L-1);    % size of the matrices
   h = 2^(L-1) / N;    % size of the interval

   residualout = 0*v; % zero in first and last element for bc
   residualout(2:n) = f(2:n) + (v(1:n-1) + v(3:n+1) - 2*v(2:n) )/h^2;
end

