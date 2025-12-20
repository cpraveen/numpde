T  = 1.0;     % final time
dt = 0.01;    % initial time step
N  = 50;      % initial grid points
Mx = N
My = N
ngrid = 4;   % no of grids for conv study

err = zeros(ngrid,1);
for i=1:ngrid  % loop over different grids
   err(i) = peaceman_rachford(T, Mx, My, dt, false);
   Mx = 2*Mx;
   My = 2*My;
   dt = 0.5*dt;
end

% Compute conv rate
fprintf(1,'%d  %e\n', 1, err(1));
for i=2:ngrid
   rate = log(err(i-1)/err(i)) / log(2);
   fprintf(1,'%d  %e %f\n', i, err(i), rate);
end
