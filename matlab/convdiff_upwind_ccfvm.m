% convection-diffusion using upwind cell-centered fvm
% -(1/Pe) u'' + u' = 0 in (0,1)
% u(0) = a   and  u(1) = b
function convdiff_upwind_ccfvm(Pe, N)
a = 0; b = 1; h = 1/N; P = h*Pe;
A = zeros(N,N);
A(1,1) = 1+2/P; A(1,2) = -1/P;
for i=2:N-1
   A(i,i-1) = -(1+1/P);
   A(i,i)   = 1+2/P;
   A(i,i+1) = -1/P;
end
A(N,N-1) = -(1+1/P); A(N,N) = 1 + 3/P;
F=zeros(N,1);
F(1) = (1+1/P)*a; F(N) = 2/P*b;
u = A \ F;
x = linspace(0.5*h, 1-0.5*h, N);
xe = linspace(0,1,100);
ue = a + (b-a)*(exp((xe-1)*Pe) - exp(-Pe))/(1 - exp(-Pe));
plot(xe,ue,'-',x,u,'o--','LineWidth',2)
legend('Exact', 'Upwind CCFVM')
