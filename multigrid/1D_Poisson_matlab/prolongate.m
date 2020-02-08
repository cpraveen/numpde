%------------------------------------------------------
% prolongate.m

function prolonged = prolongate(v, L)
   global N;
   n = N / 2^(L-1);    % size of the matrices
   c=zeros(1,2*n +1);
   c(1:2:end) = v;     % entries of v are transferred as they are
   c(2:2:end) = ( v(1:n) + v(2:n+1) )/2;
   prolonged = c;
end
