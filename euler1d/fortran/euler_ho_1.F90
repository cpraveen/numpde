! FVM for Euler 1d
!    flux   = lxf
!    recon  = first, minmod, wenojs, wenoz
!    test   =  use include file
module constants
   implicit none
   integer,parameter :: nvar = 3
   integer,parameter :: ilxf = 1
   integer,parameter :: ifirst = 1, iminmod = 2, iwenojs = 3, iwenoz = 4
   integer,parameter :: iNeumann = 1,  iWall = 2
   integer :: iflux, irecon
   integer :: ncel, nbeg, nend
   real    :: cfl
end module constants

! Include file to describe test case
include 'sod.F90'

! Convert primitive variables to conserved variables
subroutine prim2con(v,u)
   use constants
   use TestData, only : gam
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
   use TestData, only : gam
   implicit none
   real,intent(in)    :: u(nvar)
   real,intent(inout) :: v(nvar)

   v(1) = u(1)
   v(2) = u(2) / u(1)
   v(3) = (gam - 1.0) * (u(3) - 0.5 * u(2)**2 / u(1))
end subroutine con2prim

! Maximum absolute speed for euler = |vel| + sound_speed
subroutine max_speed(v, s)
   use constants
   use TestData, only : gam
   implicit none
   real,intent(in) :: v(nvar)
   real,intent(inout) :: s

   s = abs(v(2)) + sqrt(gam*v(3)/v(1))
end subroutine max_speed

! Compute euler flux from primitive variables v
subroutine euler_flux(v, flux)
   use constants
   use TestData, only : gam
   implicit none
   real,intent(in)    :: v(nvar)
   real,intent(inout) :: flux(nvar)
   ! Local variables
   real :: E

   flux(1) = v(1) * v(2)
   flux(2) = v(3) + v(1) * v(2)**2
   E       = v(3) / (gam - 1.0) + 0.5 * v(1) * v(2)**2
   flux(3) = (E + v(3))*v(2)
end subroutine euler_flux

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
   real,parameter :: beta = 1.0  ! beta in [1,2]
   real,external  :: minmod

   do c=1,n
      db = uj(c) - ujm1(c)
      dc = 0.5*(ujp1(c) - ujm1(c))
      df = ujp1(c) - uj(c)
      u(c) = uj(c) + 0.5 * minmod(beta*db, dc, beta*df)
   enddo

end subroutine reconstruct_minmod

! Classical WENO of Jiang and Shu
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

! Improved WENO at extrema
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

   if(irecon == ifirst)then
      u = uj
   else if(irecon == iminmod)then
      call reconstruct_minmod(nvar, ujm1, uj, ujp1, u)
   else if(irecon == iwenojs)then
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
   implicit none
   real,intent(in)    :: u(nvar,nbeg:nend)
   real,intent(inout) :: res(nvar,0:ncel+1)
   ! Local variables
   integer :: i
   real    :: ul(nvar), ur(nvar), flux(nvar)

   res = 0.0

   ! All faces
   do i=0,ncel
      ! Face is between i and i+1
      call reconstruct(u(:,i-2), u(:,i-1), u(:,i), u(:,i+1), u(:,i+2), ul)
      call reconstruct(u(:,i+3), u(:,i+2), u(:,i+1), u(:,i), u(:,i-1), ur)
      call num_flux(ul, ur, flux)
      res(:,i  ) = res(:,i  ) + flux
      res(:,i+1) = res(:,i+1) - flux
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

! Read parameters from file input.txt
subroutine read_input()
   use constants
   integer :: fid = 10
   character(len=20) :: recon

   open(fid,file='input.txt',status='old')
   read(fid,*) cfl
   read(fid,*) ncel
   read(fid,*) recon
   close(fid)

   if(recon == 'first')then
      irecon = ifirst
   else if(recon == 'minmod')then
      irecon = iminmod
   else if(recon == 'wenojs')then
      irecon = iwenojs
   else if(recon == 'wenoz')then
      irecon = iwenoz
   else
      stop 'Unknown recon in input.txt'
   endif

   write(*,*) 'cfl      = ', cfl
   write(*,*) 'no cells = ', ncel
   write(*,*) 'recon    = ', trim(recon)

end subroutine read_input

! This where execution starts
program main
   use constants
   use TestData
   implicit none
   real,allocatable  :: u0(:,:), u(:,:), res(:,:), xc(:)
   real              :: dx, dt, t
   integer           :: i, iter
   real,external     :: compute_dt

   iflux = ilxf
   call read_input()

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

      ! Stage 1 of SSPRK3
      call compute_residual(u, res)
      u(:,1:ncel)    = u(:,1:ncel) - (dt/dx)*res(:,1:ncel)
      call fill_ghost(u)

      ! Stage 2 of SSPRK3
      call compute_residual(u, res)
      u(:,1:ncel)    = (3.0/4.0)*u0(:,1:ncel) + &
                       (1.0/4.0)*(u(:,1:ncel) - (dt/dx)*res(:,1:ncel))
      call fill_ghost(u)

      ! Stage 3 of SSPRK3
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
