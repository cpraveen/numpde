module constants
   implicit none
   integer,parameter :: nvar = 4
   real,parameter :: gam = 1.4
   integer,parameter :: ifirst = 1, iminmod = 2, iwenojs = 3, iwenoz = 4
   integer :: irecon
   integer :: nx, ny, nx1, nx2, ny1, ny2
   real    :: xmin, xmax, ymin, ymax, dx, dy
end module constants

module WorkData
   implicit none
   real,allocatable :: v(:,:,:), fx(:,:,:), fy(:,:,:)
end module Workdata

! Convert primitive variables to conserved variables
subroutine prim2con(v,u)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: u(nvar)

   u(1) = v(1)
   u(2) = v(1)*v(2)
   u(3) = v(1)*v(3)
   u(4) = v(4)/(gam-1.0) + 0.5*v(1)*(v(2)**2 + v(3)**2)
end subroutine prim2con

! Convert primitive variables to conserved variables
subroutine con2prim(u, v)
   use constants
   implicit none
   real,intent(in)    :: u(nvar)
   real,intent(inout) :: v(nvar)

   v(1) = u(1)
   v(2) = u(2) / u(1)
   v(3) = u(3) / u(1)
   v(4) = (gam-1.0) * (u(4) - 0.5 * (u(2)**2 + u(3)**2) / u(1))
end subroutine con2prim

! Maximum absolute speed for euler = |vel| + sound_speed
subroutine max_speed_x(v, s)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: s

   s = abs(v(2)) + sqrt(gam*v(4)/v(1))
end subroutine max_speed_x

subroutine max_speed_y(v, s)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: s

   s = abs(v(3)) + sqrt(gam*v(4)/v(1))
end subroutine max_speed_y

subroutine euler_fluxes(v, fx, fy)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: fx(nvar), fy(nvar)
   ! Local variables
   real :: E

   E     = v(4)/(gam-1.0) + 0.5*v(1)*(v(2)**2 + v(3)**2)

   fx(1) = v(1)*v(2)
   fx(2) = v(4) + v(1)*v(2)**2
   fx(3) = v(1)*v(2)*v(3)
   fx(4) = (E + v(4))*v(2)

   fy(1) = v(1)*v(3)
   fy(2) = v(1)*v(2)*v(3)
   fy(3) = v(4) + v(1)*v(3)**2
   fy(4) = (E + v(4))*v(3)
end subroutine euler_fluxes

! Initial condition as a function of x
subroutine initial_condition(x, y, u)
   use constants
   implicit none
   real,intent(in)    :: x, y
   real,intent(inout) :: u(nvar)
   ! Local variables
   real :: M, alpha, beta, r2, v(nvar)
   real,parameter :: M_PI = 4.0*atan(1.0)

   ! Set initial condition for isentropic vortex
   M = 0.5
   alpha = 45.0
   beta = 5.0
   r2 = x*x + y*y
   v(1) =  (1.0 - (gam-1.0)*(beta*beta)/(8.0*gam*M_PI*M_PI)*exp(1.0-r2))**(1.0/(gam-1.0))
   v(2) =  M*cos(alpha*M_PI/180.0) - beta/(2.0*M_PI)*y*exp(0.5*(1.0-r2))
   v(3) =  M*sin(alpha*M_PI/180.0) + beta/(2.0*M_PI)*x*exp(0.5*(1.0-r2))
   v(4) =  v(1)**gam

   call prim2con(v, u)
end subroutine initial_condition

! Compute time step from cfl condition
real function compute_dt(u)
   use constants
   implicit none
   real,intent(in) :: u(nvar,nx1:nx2,ny1:ny2)
   ! Local variables
   integer :: i, j
   real    :: dt, v(nvar), sx, sy

   dt = 1.0e20
   do j=1,ny
      do i=1,nx
         call con2prim(u(:,i,j), v)
         call max_speed_x(v, sx)
         call max_speed_y(v, sy)
         dt = min(dt, 1.0/(sx/dx + sy/dy))
      enddo
   enddo

   compute_dt = dt

end function compute_dt

! Fill 3 ghost cells on each side using periodic bc
subroutine fill_ghost(u)
   use constants
   implicit none
   real,intent(inout) :: u(nvar,nx1:nx2,ny1:ny2)
   ! Local variables
   integer :: i, j

   do j=1,ny
      ! Ghost cells on left side
      u(:, 0,j) = u(:,nx  ,j)
      u(:,-1,j) = u(:,nx-1,j)
      u(:,-2,j) = u(:,nx-2,j)

      ! Ghost cells on right side
      u(:,nx+1,j) = u(:,1,j)
      u(:,nx+2,j) = u(:,2,j)
      u(:,nx+3,j) = u(:,3,j)
   enddo

   do i=1,nx
      ! ghost cells on bottom side
      u(:,i, 0) = u(:,i,ny  )
      u(:,i,-1) = u(:,i,ny-1)
      u(:,i,-2) = u(:,i,ny-2)

      ! Ghost cells on top side
      u(:,i,ny+1) = u(:,i,1)
      u(:,i,ny+2) = u(:,i,2)
      u(:,i,ny+3) = u(:,i,3)
   enddo

end subroutine fill_ghost

subroutine reconstruct_wenojs(n, ujm2, ujm1, uj, ujp1, ujp2, u)
   implicit none
   integer,intent(in) :: n
   real,intent(in)    :: ujm2(n), ujm1(n), uj(n), ujp1(n), ujp2(n)
   real,intent(inout) :: u(n)
   ! Local variables
   real,parameter :: eps = 1.0e-6
   real,parameter :: gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0
   real,dimension(n) :: beta1, beta2, beta3, w1, w2, w3, u1, u2, u3

   beta1 = (13.0/12.0)*(ujm2 - 2.0*ujm1 + uj)**2 + &
           (1.0/4.0)*(ujm2 - 4.0*ujm1 + 3.0*uj)**2
   beta2 = (13.0/12.0)*(ujm1 - 2.0*uj + ujp1)**2 + &
           (1.0/4.0)*(ujm1 - ujp1)**2
   beta3 = (13.0/12.0)*(uj - 2.0*ujp1 + ujp2)**2 + &
           (1.0/4.0)*(3.0*uj - 4.0*ujp1 + ujp2)**2

   w1 = gamma1 / (eps+beta1)**2
   w2 = gamma2 / (eps+beta2)**2
   w3 = gamma3 / (eps+beta3)**2

   u1 = (1.0/3.0)*ujm2 - (7.0/6.0)*ujm1 + (11.0/6.0)*uj
   u2 = -(1.0/6.0)*ujm1 + (5.0/6.0)*uj + (1.0/3.0)*ujp1
   u3 = (1.0/3.0)*uj + (5.0/6.0)*ujp1 - (1.0/6.0)*ujp2

   u = (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

end subroutine reconstruct_wenojs

subroutine reconstruct_wenoz(n, ujm2, ujm1, uj, ujp1, ujp2, u)
   implicit none
   integer,intent(in) :: n
   real,intent(in)    :: ujm2(n), ujm1(n), uj(n), ujp1(n), ujp2(n)
   real,intent(inout) :: u(n)
   ! Local variables
   real,parameter :: eps = 1.0e-6
   real,parameter :: gamma1=1.0/10.0, gamma2=3.0/5.0, gamma3=3.0/10.0
   real,dimension(n) :: beta1, beta2, beta3, tau, w1, w2, w3, u1, u2, u3

   beta1 = (13.0/12.0)*(ujm2 - 2.0*ujm1 + uj)**2 + &
           (1.0/4.0)*(ujm2 - 4.0*ujm1 + 3.0*uj)**2
   beta2 = (13.0/12.0)*(ujm1 - 2.0*uj + ujp1)**2 + &
           (1.0/4.0)*(ujm1 - ujp1)**2
   beta3 = (13.0/12.0)*(uj - 2.0*ujp1 + ujp2)**2 + &
           (1.0/4.0)*(3.0*uj - 4.0*ujp1 + ujp2)**2

   tau = abs(beta1 - beta3)
   w1 = gamma1 * (1.0 + (tau / (eps+beta1))**2)
   w2 = gamma2 * (1.0 + (tau / (eps+beta2))**2)
   w3 = gamma3 * (1.0 + (tau / (eps+beta3))**2)

   u1 = (1.0/3.0)*ujm2 - (7.0/6.0)*ujm1 + (11.0/6.0)*uj
   u2 = -(1.0/6.0)*ujm1 + (5.0/6.0)*uj + (1.0/3.0)*ujp1
   u3 = (1.0/3.0)*uj + (5.0/6.0)*ujp1 - (1.0/6.0)*ujp2

   u = (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

end subroutine reconstruct_wenoz

! Reconstruct left state at interface between j and j+1
subroutine reconstruct(ujm2, ujm1, uj, ujp1, ujp2, u)
   use constants
   implicit none
   real,dimension(nvar),intent(in) :: ujm2, ujm1, uj, ujp1, ujp2
   real,dimension(nvar),intent(inout) :: u

   if(irecon == iwenojs)then
      call reconstruct_wenojs(nvar, ujm2, ujm1, uj, ujp1, ujp2, u)
   else if(irecon == iwenoz)then
      call reconstruct_wenoz(nvar, ujm2, ujm1, uj, ujp1, ujp2, u)
   else
      stop 'Unknown value of irecon'
   endif

end subroutine reconstruct

! Compute finite volume residual R
! dx*du/dt + R = 0
subroutine compute_residual(u, res)
   use constants
   use WorkData
   implicit none
   real,intent(in)    :: u(nvar,nx1:nx2,ny1:ny2)
   real,intent(inout) :: res(nvar,0:nx+1,0:ny+1)
   ! Local variables
   integer :: i, j
   real    :: fim2(nvar), fim1(nvar), fi(nvar), fip1(nvar), fip2(nvar), fip3(nvar)
   real    :: lamx1, lamx2, lamx, lamy1, lamy2, lamy
   real    :: fm(nvar), fp(nvar), flux(nvar)

   do j=ny1,ny2
      do i=nx1,nx2
         call con2prim(u(:,i,j), v(:,i,j))
         call euler_fluxes(v(:,i,j), fx(:,i,j), fy(:,i,j))
      enddo
   enddo

   res = 0.0

   ! Vertical faces
   do j=1,ny
      do i=0,nx
         ! Face is between i and i+1
         call max_speed_x(v(:,i  ,j), lamx1)
         call max_speed_x(v(:,i+1,j), lamx2)
         lamx = max(lamx1, lamx2)

         ! Positive flux
         fim2 = 0.5*(fx(:,i-2,j) + lamx * u(:,i-2,j))
         fim1 = 0.5*(fx(:,i-1,j) + lamx * u(:,i-1,j))
         fi   = 0.5*(fx(:,i  ,j) + lamx * u(:,i  ,j))
         fip1 = 0.5*(fx(:,i+1,j) + lamx * u(:,i+1,j))
         fip2 = 0.5*(fx(:,i+2,j) + lamx * u(:,i+2,j))

         call reconstruct(fim2, fim1, fi, fip1, fip2, fp)

         ! Negative flux
         fim1 = 0.5*(fx(:,i-1,j) - lamx * u(:,i-1,j))
         fi   = 0.5*(fx(:,i  ,j) - lamx * u(:,i  ,j))
         fip1 = 0.5*(fx(:,i+1,j) - lamx * u(:,i+1,j))
         fip2 = 0.5*(fx(:,i+2,j) - lamx * u(:,i+2,j))
         fip3 = 0.5*(fx(:,i+3,j) - lamx * u(:,i+3,j))

         call reconstruct(fip3, fip2, fip1, fi, fim1, fm)

         flux = (fp + fm)*dy
         res(:,i  ,j) = res(:,i  ,j) + flux
         res(:,i+1,j) = res(:,i+1,j) - flux
      enddo
   enddo

   ! Horizontal faces
   do j=0,ny
      do i=1,nx
         ! Face is between j and j+1
         call max_speed_y(v(:,i,j  ), lamy1)
         call max_speed_y(v(:,i,j+1), lamy2)
         lamy = max(lamy1, lamy2)

         ! Positive flux
         fim2 = 0.5*(fy(:,i,j-2) + lamy * u(:,i,j-2))
         fim1 = 0.5*(fy(:,i,j-1) + lamy * u(:,i,j-1))
         fi   = 0.5*(fy(:,i,j  ) + lamy * u(:,i,j  ))
         fip1 = 0.5*(fy(:,i,j+1) + lamy * u(:,i,j+1))
         fip2 = 0.5*(fy(:,i,j+2) + lamy * u(:,i,j+2))

         call reconstruct(fim2, fim1, fi, fip1, fip2, fp)

         ! Negative flux
         fim1 = 0.5*(fy(:,i,j-1) - lamy * u(:,i,j-1))
         fi   = 0.5*(fy(:,i,j  ) - lamy * u(:,i,j  ))
         fip1 = 0.5*(fy(:,i,j+1) - lamy * u(:,i,j+1))
         fip2 = 0.5*(fy(:,i,j+2) - lamy * u(:,i,j+2))
         fip3 = 0.5*(fy(:,i,j+3) - lamy * u(:,i,j+3))

         call reconstruct(fip3, fip2, fip1, fi, fim1, fm)

         flux = (fp + fm)*dx
         res(:,i,j  ) = res(:,i,j  ) + flux
         res(:,i,j+1) = res(:,i,j+1) - flux
      enddo
   enddo

end subroutine compute_residual

subroutine getfilename(filename, ext)
   implicit none
   character(len=*) :: filename, ext
   ! Local variables
   integer,save :: fileid = 0

   if(fileid <= 9)then
      write(unit=filename, fmt='(A5,I1)') trim(filename)//'00', fileid
   else if(fileid <= 99)then
      write(unit=filename, fmt='(A4,I2)') trim(filename)//'0', fileid
   else if(fileid <= 999)then
      write(unit=filename, fmt='(A3,I3)') trim(filename), fileid
   else
      stop 'getfilename: counter is too large'
   endif

   filename = trim(filename)//'.'//trim(ext)

   fileid = fileid + 1

end subroutine getfilename

! Save solution to file
subroutine savesol(time, u)
   use constants
   implicit none
   real,intent(in)    :: time, u(nvar,nx1:nx2,ny1:ny2)
   ! Local variables
   integer :: i, j, fid
   real    :: x, y, v(nvar)
   character(len=512) :: filename

   fid = 10
   filename = 'sol'
   call getfilename(filename, 'plt')
   open(fid,file=trim(filename))
   write(fid,*)'TITLE = "Euler"'
   write(fid,*)'VARIABLES = "x","y","rho","vx","vy","Pre"'
   write(fid,*)'ZONE STRANDID=1, SOLUTIONTIME=',time,', I=',nx,', J=',ny,&
              ', DATAPACKING=POINT'
   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         call con2prim(u(:,i,j), v)
         write(fid,*) x, y, v(1), v(2), v(3), v(4)
      enddo
   enddo
   close(fid)
   print*,'Saved solution into ',trim(filename)

end subroutine savesol

program main
   use constants
   use WorkData
   implicit none
   real,allocatable  :: u0(:,:,:), u(:,:,:), res(:,:,:)
   real              :: cfl, dt, t, Tf, x, y
   integer           :: i, j, iter
   real,external     :: compute_dt

   xmin = -5.0
   xmax =  5.0
   ymin = -5.0
   ymax =  5.0
   Tf   = 10.0 * sqrt(2.0) / 0.5

   cfl    = 0.95
   nx     = 100
   ny     = 100
   irecon = iwenoz

   ! 3 ghost cells on each side, needed for WENO
   nx1 = -2
   nx2 = nx + 3
   ny1 = -2
   ny2 = ny + 3
   print*,'nx,ny =', nx, ny

   ! Allocate memory
   allocate(u(nvar,nx1:nx2,ny1:ny2), &
            u0(nvar,nx1:nx2,ny1:ny2), &
            res(nvar,0:nx+1,0:ny+1))
   allocate(v(nvar,nx1:nx2,ny1:ny2), &
            fx(nvar,nx1:nx2,ny1:ny2), &
            fy(nvar,nx1:nx2,ny1:ny2))

   ! Cell size
   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny
   print*,'dx, dy = ', dx, dy

   ! Set initial condition
   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         call initial_condition(x, y, u(:,i,j))
      enddo
   enddo
   call fill_ghost(u)

   t    = 0.0
   iter = 0
   call savesol(t, u)
   do while(t < Tf)
      u0 = u
      dt = cfl * compute_dt(u)
      if(t+dt > Tf) dt = Tf - t

      ! Stage 1
      call compute_residual(u, res)
      u(:,1:nx,1:ny)    = u(:,1:nx,1:ny) - (dt/(dx*dy))*res(:,1:nx,1:ny)
      call fill_ghost(u)

      ! Stage 2
      call compute_residual(u, res)
      u(:,1:nx,1:ny) = (3.0/4.0)*u0(:,1:nx,1:ny) + &
                       (1.0/4.0)*(u(:,1:nx,1:ny) - (dt/(dx*dy))*res(:,1:nx,1:ny))
      call fill_ghost(u)

      ! Stage 3
      call compute_residual(u, res)
      u(:,1:nx,1:ny) = (1.0/3.0)*u0(:,1:nx,1:ny) + &
                       (2.0/3.0)*(u(:,1:nx,1:ny) - (dt/(dx*dy))*res(:,1:nx,1:ny))
      call fill_ghost(u)

      t    = t + dt
      iter = iter + 1
      print*,iter,dt,t

      if(mod(iter,100) == 0) call savesol(t, u)
   enddo

end program main