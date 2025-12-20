% Construct Laplace FD matrix in 2D
n=6;   % grid points in each direction
m=n-2; % interior grid points in each direction
I = speye(m); 
e = ones(m,1);
D = spdiags([e -2*e e],[-1 0 1],m,m);
A = -(kron(I,D) + kron(D,I));
full(A)
spy(A,'s'); axis('equal')
