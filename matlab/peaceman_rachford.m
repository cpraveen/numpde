% Solve 2-d heat equation using Peaceman-Rachford scheme
% ut = uxx + uyy + F
% Exact solution: u = 1/[1 + 50(x-x0)^2 + 50(y-y0)^2]
%   x0(t) = 1/2 + 1/4*cos(2*pi*t)
%   y0(t) = 1/2 + 1/4*sin(2*pi*t)
% Input parameters:
%   Mx = No. of grid points along x
%   My = No. of grid points along y
%    T = Final time
function error = peaceman_rachford(Mx, My, T)

mu = 1;
dx = 1/(Mx-1);
dy = 1/(My-1);

dt = 0.01;

rx = mu*dt/dx^2; ry = mu*dt/dy^2;
fprintf(1,'rx,ry = %e, %e\n', rx, ry);
ex = ones(Mx-2,1);
Ax = spdiags([-0.5*rx*ex, (1+rx)*ex, -0.5*rx*ex], -1:1, Mx-2, Mx-2);
ey = ones(My-2,1);
Ay = spdiags([-0.5*ry*ey, (1+ry)*ey, -0.5*ry*ey], -1:1, My-2, My-2);

x  = linspace(0,1,Mx);
y  = linspace(0,1,My);
[X,Y] = ndgrid(x,y);

u    = zeros(Mx,My);
unph = zeros(Mx,My);

u = ExactSolution(X,Y,0);

contour(X,Y,u), axis square, pause(2)

t = 0;
while t < T
   if t+dt > T, dt = T - t; end
   % Set boundary conditions
   unph(1:Mx,1   ) = ExactSolution(X(1:Mx,1   ),Y(1:Mx,1   ),t+0.5*dt);
   unph(1:Mx,My  ) = ExactSolution(X(1:Mx,My  ),Y(1:Mx,My  ),t+0.5*dt);
   unph(1,   1:My) = ExactSolution(X(1,   1:My),Y(1,   1:My),t+0.5*dt);
   unph(Mx,  1:My) = ExactSolution(X(Mx,  1:My),Y(Mx,  1:My),t+0.5*dt);
   for k=2:My-1
      rhs = 0.5*dt*SourceTerm(X(2:Mx-1,k), Y(2:Mx-1,k),t+0.5*dt);
      rhs = rhs + u(2:Mx-1,k) + ...
            0.5*ry*(u(2:Mx-1,k-1) - 2*u(2:Mx-1,k) + u(2:Mx-1,k+1));
      rhs(1) = rhs(1) + 0.5*rx*unph(1,k);
      rhs(end) = rhs(end) + 0.5*rx*unph(Mx,k);
      unph(2:Mx-1,k) = Ax \ rhs;
   end
   % Set boundary conditions
   u(1:Mx,1   ) = ExactSolution(X(1:Mx,1   ),Y(1:Mx,1   ),t+dt);
   u(1:Mx,My  ) = ExactSolution(X(1:Mx,My  ),Y(1:Mx,My  ),t+dt);
   u(1,   1:My) = ExactSolution(X(1,   1:My),Y(1,   1:My),t+dt);
   u(Mx,  1:My) = ExactSolution(X(Mx,  1:My),Y(Mx,  1:My),t+dt);
   for j=2:Mx-1
      rhs = 0.5*dt*SourceTerm(X(j,2:My-1), Y(j,2:My-1),t+0.5*dt);
      rhs = rhs + unph(j,2:My-1) + ...
            0.5*rx*(unph(j-1,2:My-1) - 2*unph(j,2:My-1) + unph(j+1,2:My-1));
      rhs(1) = rhs(1) + 0.5*ry*u(j,1);
      rhs(end) = rhs(end) + 0.5*ry*u(j,My);
      u(j,2:My-1) = Ay \ rhs';
   end
   umin=min(min(u)); umax=max(max(u));
   t = t + dt; fprintf(1,'t=%e, min,max=%e, %e\n',t,umin,umax);
   st=strcat('t = ', num2str(t));
   contour(X,Y,u), axis square, title(st), pause(.1)
end

% Compute L2 norm of error
ue = ExactSolution(X,Y,T);
error = norm(u-ue,'fro');
error = error*sqrt(dx*dy);

end
%-------------------------------------------------------------------------------
% Exact solution
%-------------------------------------------------------------------------------
function u = ExactSolution(x,y,t)

x0 = 0.5 + 0.25 * cos(2*pi*t);
y0 = 0.5 + 0.25 * sin(2*pi*t);
u  = 1./(1 + 50*(x-x0).^2 + 50*(y-y0).^2);

end
%-------------------------------------------------------------------------------
% Source term
%-------------------------------------------------------------------------------
function s = SourceTerm(x,y,t)

x0 = 0.5 + 0.25 * cos(2*pi*t);
y0 = 0.5 + 0.25 * sin(2*pi*t);
u  = 1./(1 + 50*(x-x0).^2 + 50*(y-y0).^2);

x0t = -0.25 * sin(2*pi*t) * (2*pi);
y0t =  0.25 * cos(2*pi*t) * (2*pi);
 
ut = -u.^2 .* (-100*(x-x0)*x0t - 100*(y-y0)*y0t);

ux = -u.^2 .* 100.*(x-x0);
uy = -u.^2 .* 100.*(y-y0);

uxx = -u.^2 * 100 - 200*u.*ux.*(x-x0);
uyy = -u.^2 * 100 - 200*u.*uy.*(y-y0);

s   = ut - uxx - uyy;
end
