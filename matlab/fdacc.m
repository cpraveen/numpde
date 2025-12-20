% f(x) = sin(x), at x = pi/4,
% Compute error in f'(x) using forward difference and central difference
clear all
f = @(x) sin(x);
x = pi/4;
df= cos(x); % exact derivative
dx=0.01;
res=[];
for j=1:40
   h(j) = dx;
   fd = ( f(x+h(j)) - f(x) ) / h(j);
   efd(j) = abs(fd - df);
   cd = ( f(x+h(j)) - f(x-h(j)) ) / (2*h(j));
   ecd(j) = abs(cd - df);
   if j>1
      p1 = log( efd(j-1)/efd(j) ) / log(2);
      p2 = log( ecd(j-1)/ecd(j) ) / log(2);
      res=[res; h(j), efd(j), p1, ecd(j), p2];
   end
   dx = dx/2;
end
format shorte; res
loglog(h, efd, 'o-', h, ecd, '*-')
set(gca,'FontSize', 16)
xlabel('h')
ylabel('abs(Error)')
legend('Forward difference', 'Central difference', 'Location', 'SouthWest')
print -dpdf fdacc.pdf
