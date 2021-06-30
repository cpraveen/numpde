! FVM for Euler 1d
!    flux   = lxf
!    recon  = first order
!    test   = sod
module constants
   implicit none
   integer,parameter :: nvar = 3
   real,parameter :: gam = 1.4
   integer,parameter :: ilxf = 1
   integer :: iflux
end module constants

! Convert primitive variables to conserved variables
subroutine prim2con(v,u)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: u(nvar)

   u(1) = v(1)
   u(2) = v(1) * v(2)
   u(3) = v(3) / (gam - 1.0) + 0.5 * v(1) * v(2)**2
end subroutine prim2con

! Convert primitive variables to conserved variables
subroutine con2prim(u, v)
   use constants
   implicit none
   real,intent(in)    :: u(nvar)
   real,intent(inout) :: v(nvar)

   v(1) = u(1)
   v(2) = u(2) / u(1)
   v(3) = (gam - 1.0) * (u(3) - 0.5 * u(2)**2 / u(1))
end subroutine con2prim

! Largest eigenvalue in magnitude
subroutine max_speed(v, s)
   use constants
   implicit none
   real,intent(in) :: v(nvar)
   real,intent(inout) :: s

   s = abs(v(2)) + sqrt(gam*v(3)/v(1))
end subroutine max_speed

subroutine euler_flux(v, flux)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: flux(nvar)
   ! Local variables
   real :: E

   flux(1) = v(1) * v(2)
   flux(2) = v(3) + v(1) * v(2)**2
   E       = v(3) / (gam - 1.0) + 0.5 * v(1) * v(2)**2
   flux(3) = (E + v(3)) * v(2)
end subroutine euler_flux

! Sod shock tube case
subroutine initial_condition(x, u)
   use constants
   implicit none
   real,intent(in)    :: x
   real,intent(inout) :: u(nvar)
   ! Local variables
   real :: v(nvar)

   if(x < 0.5)then
      v(1) = 1.0
      v(2) = 0.0
      v(3) = 1.0
   else
      v(1) = 0.125
      v(2) = 0.0
      v(3) = 0.1
   endif

   call prim2con(v, u)
end subroutine initial_condition

! Compute numerical flux: Rusanov flux
subroutine lxf_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), sl, sr, fl(nvar), fr(nvar), lam

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   call max_speed(vl, sl)
   call max_speed(vr, sr)
   lam = max(sl, sr)

   call euler_flux(vl, fl)
   call euler_flux(vr, fr)

   nflux = 0.5*(fl + fr) - 0.5*lam*(ur - ul)
end subroutine lxf_flux

! Compute numerical flux
subroutine num_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)

   if(iflux == ilxf)then
      call lxf_flux(ul, ur, nflux)
   else
      stop 'Only lxf_flux is available'
   endif

end subroutine num_flux

! Compute time step from cfl condition
real function compute_dt(nc, dx, u)
   use constants
   implicit none
   integer,intent(in) :: nc
   real,intent(in)    :: dx, u(nvar,nc)
   ! Local variables
   integer :: i
   real    :: dt, v(nvar), s

   dt = 1.0e20
   do i=1,nc
      call con2prim(u(:,i), v)
      call max_speed(v, s)
      dt = min(dt, dx/s)
   enddo

   compute_dt = dt

end function compute_dt

! Compute finite volume residual R
! dx*du/dt + R = 0
subroutine compute_residual(nc, u, res)
   use constants
   implicit none
   integer,intent(in) :: nc
   real,intent(in)    :: u(nvar,nc)
   real,intent(inout) :: res(nvar,nc)
   ! Local variables
   integer :: i
   real    :: flux(nvar)

   res = 0.0

   ! First face
   call num_flux(u(:,1), u(:,1), flux)
   res(:,1) = res(:,1) - flux

   ! Intermediate faces: between i and i+1
   do i=1,nc-1
      call num_flux(u(:,i), u(:,i+1), flux)
      res(:,i  ) = res(:,i  ) + flux
      res(:,i+1) = res(:,i+1) - flux
   enddo

   ! Last face
   call num_flux(u(:,nc), u(:,nc), flux)
   res(:,nc) = res(:,nc) + flux

end subroutine compute_residual

! Save solution to file
subroutine savesol(nc, xc, u)
   use constants
   implicit none
   integer,intent(in) :: nc
   real,intent(in)    :: xc(nc), u(nvar,nc)
   ! Local variables
   integer :: i, fid
   real    :: v(nvar)

   fid = 10
   open(fid, file='sol.txt')
   do i=1,nc
      call con2prim(u(:,i), v)
      write(fid,*) xc(i), v(1), v(2), v(3)
   enddo
   close(fid)
   print*,'Saved solution into sol.txt'

end subroutine savesol

! This is where program starts
program main
   use constants
   implicit none
   integer,parameter :: nc = 200
   real,parameter    :: xmin = 0.0, xmax = 1.0
   real              :: u(nvar,nc), res(nvar,nc), xc(nc)
   real              :: dx, dt, t, Tf
   integer           :: i, iter
   real,external     :: compute_dt

   dx = (xmax - xmin)/nc
   Tf = 0.2
   iflux = ilxf

   ! Make grid
   ! We locate cell faces at xmin and xmax
   do i=1,nc
      xc(i) = xmin + (i-1)*dx + 0.5*dx
   enddo

   ! Set initial condition
   do i=1,nc
      call initial_condition(xc(i), u(:,i))
   enddo

   t    = 0.0
   iter = 0
   do while(t < Tf)
      dt = compute_dt(nc, dx, u)
      if(t+dt > Tf) dt = Tf - t
      call compute_residual(nc, u, res)
      u    = u - (dt/dx)*res
      t    = t + dt
      iter = iter + 1
      print*,iter,t,dt
   enddo

   call savesol(nc, xc, u)

end program main