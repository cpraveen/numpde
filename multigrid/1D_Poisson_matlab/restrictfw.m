%------------------------------------------------------
% restrictfw.m

function restrictionout = restrictfw (r, L)
    % restriction by full weighting
    global N;
        n = N / 2^(L-1);    % size of the matrices
        c = r(1:2:end);  % c has size n/2 +1
    tempans = ( r(2:2:n-2) + r(4:2:n) + 2 * r(3:2:n-1) )/4;
    c(2:n/2) = tempans;
    restrictionout = c;
% end
