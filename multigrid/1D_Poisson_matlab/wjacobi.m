%------------------------------------------------------
% wjacobi.m

function wjreturn = wjacobi(v, f, k, L)
    global N;

    n = N / 2^(L-1);    % size of the matrices
    w = 2/3;            % weight
    h = 2^(L-1) / N;    % size of the interval

    for i = 1:k
        tempans = .5 * ( v(1:n-1) + v(3:n+1) + h*h*f(2:n) );
        vtemp = v;
        vtemp(2:n) = tempans;                   % keep boundary values unchanged
        v = (1-w)*v + w*vtemp;                      % new iterate (weighted jacobi)
    end
    wjreturn = v;
% end
