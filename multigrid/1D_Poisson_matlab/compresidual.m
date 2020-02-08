%------------------------------------------------------
% compresidual.m

function residualout = compresidual (v, f, L)
    global N;
        n = N / 2^(L-1);    % size of the matrices
        h = 2^(L-1) / N;    % size of the interval
    vtemp = v;
    tempans = (v(1:n-1) + v(3:n+1) - 2*v(2:n) )/h/h;
    vtemp(2:n) = tempans;

    residualout = f + vtemp;
% end

