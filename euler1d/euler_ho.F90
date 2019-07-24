module constants
   implicit none
   integer,parameter :: nvar = 3
   real,parameter :: gam = 1.4
   integer,parameter :: ilxf = 1
   integer,parameter :: ifirst = 1, iminmod = 2
   integer :: iflux, irecon
   integer :: ncel, nbeg, nend
end module constants

! Convert primitive variables to conserved variables
subroutine prim2con(v,u)
   use constants
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: u(nvar)

   u(1) = v(1)
   u(2) = v(1)*v(2)
   u(3) = v(3)/(gam-1.0) + 0.5*v(1)*v(2)**2
end subroutine prim2con

! Convert primitive variables to conserved variables
subroutine con2prim(u, v)
   use constants
   implicit none
   real,intent(in)    :: u(nvar)
   real,intent(inout) :: v(nvar)

   v(1) = u(1)
   v(2) = u(2) / u(1)
   v(3) = (gam-1.0) * (u(3) - 0.5 * u(2)**2 / u(1))
end subroutine con2prim

! Maximum absolute speed for euler = |vel| + sound_speed
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

   flux(1) = v(1)*v(2)
   flux(2) = v(3) + v(1)*v(2)**2
   E       = v(3)/(gam-1.0) + 0.5*v(1)*v(2)**2
   flux(3) = (E + v(3))*v(2)
end subroutine euler_flux

! Initial condition as a function of x
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

! Compute numerical flux
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
      stop 'Unknown value of iflux'
   endif

end subroutine num_flux

! Compute time step from cfl condition
real function compute_dt(dx, u)
   use constants
   implicit none
   real,intent(in)    :: dx, u(nvar,nbeg:nend)
   ! Local variables
   integer :: i
   real    :: dt, v(nvar), s

   dt = 1.0e20
   do i=1,ncel
      call con2prim(u(:,i), v)
      call max_speed(v, s)
      dt = min(dt, dx/s)
   enddo

   compute_dt = dt

end function compute_dt

! Fill 3 ghost cells on each side using neumann bc
subroutine fill_ghost(u)
   use constants
   implicit none
   real,intent(inout) :: u(nvar,nbeg:nend)

   u(:,-2) = u(:,1)
   u(:,-1) = u(:,1)
   u(:, 0) = u(:,1)

   u(:,ncel+1) = u(:,ncel)
   u(:,ncel+2) = u(:,ncel)
   u(:,ncel+3) = u(:,ncel)

end subroutine fill_ghost

! Returns minmod of three variables
real function minmod(a, b, c)
   implicit none
   real,intent(in) :: a, b, c
   ! Local variables
   real :: sa, sb, sc

   sa = sign(1.0, a)
   sb = sign(1.0, b)
   sc = sign(1.0, c)

   if(sa == sb .and. sb == sc)then
      minmod = sa * min(abs(a), abs(b), abs(c))
   else
      minmod = 0.0
   endif

end function minmod

! Minmod function
subroutine reconstruct_minmod(n, ujm1, uj, ujp1, u)
   implicit none
   integer,intent(in) :: n
   real,intent(in)    :: ujm1(n), uj(n), ujp1(n)
   real,intent(inout) :: u(n)
   ! Local variables
   integer        :: c
   real           :: db, dc, df
   real,parameter :: beta = 1.0
   real,external  :: minmod

   do c=1,n
      db = uj(c) - ujm1(c)
      dc = 0.5*(ujp1(c) - ujm1(c))
      df = ujp1(c) - uj(c)
      u(c) = uj(c) + 0.5 * minmod(beta*db, dc, beta*df)
   enddo

end subroutine reconstruct_minmod

! Reconstruct left state at interface between j and j+1
subroutine reconstruct(ujm2, ujm1, uj, ujp1, ujp2, u)
   use constants
   implicit none
   real,dimension(nvar),intent(in) :: ujm2, ujm1, uj, ujp1, ujp2
   real,dimension(nvar),intent(inout) :: u

   if(irecon == ifirst)then
      u = uj
   else if(irecon == iminmod)then
      call reconstruct_minmod(nvar, ujm1, uj, ujp1, u)
   else
      stop 'Unknown value of irecon'
   endif

end subroutine reconstruct

! Compute finite volume residual R
! dx*du/dt + R = 0
subroutine compute_residual(u, res)
   use constants
   implicit none
   real,intent(in)    :: u(nvar,nbeg:nend)
   real,intent(inout) :: res(nvar,0:ncel+1)
   ! Local variables
   integer :: i
   real    :: ul(nvar), ur(nvar), flux(nvar)

   res = 0.0

   ! All faces
   do i=1,ncel+1
      ! Face is between i-1 and i
      call reconstruct(u(:,i-3), u(:,i-2), u(:,i-1), u(:,i), u(:,i+1), ul)
      call reconstruct(u(:,i+2), u(:,i+1), u(:,i), u(:,i-1), u(:,i-2), ur)
      call num_flux(ul, ur, flux)
      res(:,i-1) = res(:,i-1) + flux
      res(:,i)   = res(:,i)   - flux
   enddo

end subroutine compute_residual

! Save solution to file
subroutine savesol(xc, u)
   use constants
   implicit none
   real,intent(in)    :: xc(ncel), u(nvar,nbeg:nend)
   ! Local variables
   integer :: i, fid
   real    :: v(nvar)

   fid = 10
   open(fid, file='sol.txt')
   do i=1,ncel
      call con2prim(u(:,i), v)
      write(fid,*) xc(i), v(1), v(2), v(3)
   enddo
   close(fid)
   print*,'Saved solution into sol.txt'

end subroutine savesol

program main
   use constants
   implicit none
   real,parameter    :: xmin = 0.0, xmax = 1.0
   real,allocatable  :: u0(:,:), u(:,:), res(:,:), xc(:)
   real              :: cfl, dx, dt, t, Tf
   integer           :: i, iter
   real,external     :: compute_dt

   cfl    = 0.95
   ncel   = 200
   Tf     = 0.2
   iflux  = ilxf
   irecon = iminmod

   ! 3 ghost cells on each side, needed for WENO
   nbeg = -2
   nend = ncel + 3

   ! Allocate memory
   allocate(xc(ncel), u(nvar,nbeg:nend), u0(nvar,nbeg:nend), &
            res(nvar,0:ncel+1))

   ! Cell size
   dx = (xmax - xmin)/ncel

   ! Make grid
   ! We locate cell faces at xmin and xmax
   do i=1,ncel
      xc(i) = xmin + (i-1)*dx + 0.5*dx
   enddo

   ! Set initial condition
   do i=1,ncel
      call initial_condition(xc(i), u(:,i))
   enddo
   call fill_ghost(u)

   t    = 0.0
   iter = 0
   do while(t < Tf)
      u0 = u
      dt = cfl * compute_dt(dx, u)
      if(t+dt > Tf) dt = Tf - t

      ! Stage 1
      call compute_residual(u, res)
      u(:,1:ncel)    = u(:,1:ncel) - (dt/dx)*res(:,1:ncel)
      call fill_ghost(u)

      ! Stage 2
      call compute_residual(u, res)
      u(:,1:ncel)    = (3.0/4.0)*u0(:,1:ncel) + &
                       (1.0/4.0)*(u(:,1:ncel) - (dt/dx)*res(:,1:ncel))
      call fill_ghost(u)

      ! Stage 3
      call compute_residual(u, res)
      u(:,1:ncel)    = (1.0/3.0)*u0(:,1:ncel) + &
                       (2.0/3.0)*(u(:,1:ncel) - (dt/dx)*res(:,1:ncel))
      call fill_ghost(u)

      t    = t + dt
      iter = iter + 1
      print*,iter,t,dt
   enddo

   call savesol(xc, u)

end program main