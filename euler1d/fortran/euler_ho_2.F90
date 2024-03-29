! FVM for Euler 1d
!    flux   = lxf, roe, hll, hllc
!    recon  = first, minmod, wenojs, wenoz
!    var    = cons and char limiting
!    test   =  use include file
module constants
   implicit none
   integer,parameter :: nvar = 3
   integer,parameter :: ilxf = 1, iroe = 2, ihll = 3, ihllc = 4, ivl = 5
   integer,parameter :: ifirst = 1, iminmod = 2, iwenojs = 3, iwenoz = 4
   integer,parameter :: iNeumann = 1,  iWall = 2
   integer :: iflux, irecon, ichar
   integer :: ncel, nbeg, nend
   real    :: cfl
   logical :: efix
end module constants

! Include one file for test case definition
#if defined(SOD)
include 'sod.F90'
#elif defined(MSOD)
include 'msod.F90'
#elif defined(SHUOSHER)
include 'shuosher.F90'
#elif defined(TORO5)
include 'toro5.F90'
#endif

include 'num_flux.F90'

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
   flux(3) = (E + v(3)) * v(2)
end subroutine euler_flux

! Compute numerical flux
subroutine num_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)

   select case(iflux)
      case(ilxf)
         call lxf_flux(ul, ur, nflux)
      case(iroe)
         call roe_flux(ul, ur, nflux)
      case(ihll)
         call hll_flux(ul, ur, nflux)
      case(ihllc)
         call hllc_flux(ul, ur, nflux)
      case(ivl)
         call vanleer_flux(ul, ur, nflux)
      case default
         stop 'Unknown value of iflux'
   end select

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
   use TestData
   implicit none
   real,intent(inout) :: u(nvar,nbeg:nend)

   if(leftbc == iWall)then
       u(1,-2) = u(1,1)
       u(1,-1) = u(1,1)
       u(1, 0) = u(1,1)

       u(2,-2) = u(2,1)
       u(2,-1) = -u(2,1)
       u(2, 0) = u(2,1)

       u(3,-2) = u(3,1)
       u(3,-1) = u(3,1)
       u(3, 0) = u(3,1)
   else ! Neumann bc
       u(:,-2) = u(:,1)
       u(:,-1) = u(:,1)
       u(:, 0) = u(:,1)
   endif

   if(rightbc == iWall)then
       u(1,ncel+1) = u(1,ncel)
       u(1,ncel+2) = u(1,ncel)
       u(1,ncel+3) = u(1,ncel)

       u(2,ncel+1) = u(2,ncel)
       u(2,ncel+2) = -u(2,ncel)
       u(2,ncel+3) = u(2,ncel)

       u(3,ncel+1) = u(3,ncel)
       u(3,ncel+2) = u(3,ncel)
       u(3,ncel+3) = u(3,ncel)
   else ! Neumann bc
       u(:,ncel+1) = u(:,ncel)
       u(:,ncel+2) = u(:,ncel)
       u(:,ncel+3) = u(:,ncel)
    endif

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
   real,dimension(nvar) :: wjm2, wjm1, wj, wjp1, wjp2
   real :: R(nvar,nvar), L(nvar,nvar), v(nvar)

   ! For first order, skip everything else
   if(irecon == ifirst)then
      u = uj
      return
   endif

   if(ichar == 1)then
      ! Convert to characteristic variables
      call Eig_Vec(uj, ujp1, L, R)
      wjm2 = matmul(L, ujm2)
      wjm1 = matmul(L, ujm1)
      wj   = matmul(L, uj)
      wjp1 = matmul(L, ujp1)
      wjp2 = matmul(L, ujp2)
   else
      wjm2 = ujm2
      wjm1 = ujm1
      wj   = uj
      wjp1 = ujp1
      wjp2 = ujp2
   endif

   if(irecon == iminmod)then
      call reconstruct_minmod(nvar, wjm1, wj, wjp1, v)
   else if(irecon == iwenojs)then
      call reconstruct_wenojs(nvar, wjm2, wjm1, wj, wjp1, wjp2, v)
   else if(irecon == iwenoz)then
      call reconstruct_wenoz(nvar, wjm2, wjm1, wj, wjp1, wjp2, v)
   else
      stop 'Unknown value of irecon'
   endif

   if(ichar == 1)then
      ! Convert back to conserved variables
      u = matmul(R, v)
   else
      u = v
   endif

end subroutine reconstruct

! Form matrix of left and right eigenvectors of flux jacobian
! ul, ur = conserved variables
subroutine Eig_Vec(ul, ur, L, R)
   use constants
   use TestData
   implicit none
   real :: ul(nvar), ur(nvar), u(nvar)
   real :: R(nvar,nvar), L(nvar,nvar)
   real :: H, v(nvar), a, M

   u = 0.5*(ul + ur)
   call con2prim(u, v)
   H = (u(3) + v(3)) / u(1)
   a = sqrt(gam * v(3) / v(1))

   R(1,1) = 1.0;      R(1,2) = 1.0;         R(1,3) = 1.0
   R(2,1) = v(2)-a;   R(2,2) = v(2);        R(2,3) = v(2)+a
   R(3,1) = H-v(2)*a; R(3,2) = 0.5*v(2)**2; R(3,3) = H+v(2)*a

   M      = v(2)/a

   L(1,1) = 0.25*(gam-1.0)*M**2 + 0.5*M
   L(1,2) = -0.5*(gam-1.0)*M/a - 0.5/a
   L(1,3) = 0.5*(gam-1.0)/a**2

   L(2,1) = 1.0 - 0.5*(gam-1.0)*M**2
   L(2,2) = (gam-1.0)*M/a
   L(2,3) = -(gam-1.0)/a**2

   L(3,1) = 0.25*(gam-1.0)*M**2 - 0.5*M
   L(3,2) = -0.5*(gam-1.0)*M/a + 0.5/a
   L(3,3) = 0.5*(gam-1.0)/a**2

end subroutine Eig_Vec

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
   character(len=20) :: flux, recon, var

   open(fid,file='input.txt',status='old')
   read(fid,*) cfl
   read(fid,*) ncel
   read(fid,*) recon
   read(fid,*) var
   read(fid,*) flux
   read(fid,*) efix
   close(fid)

   if(flux == 'lxf')then
      iflux = ilxf
   else if(flux == 'roe')then
      iflux = iroe
   else if(flux == 'hll')then
      iflux = ihll
   else if(flux == 'hllc')then
      iflux = ihllc
   else if(flux == 'vl')then
      iflux = ivl
   else
      stop 'Unknown flux in input.txt'
   endif

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

   if(var == 'cons')then
      ichar = 0
   else if(var == 'char')then
      ichar = 1
   else
      stop 'Unknown var in input.txt'
   endif

   write(*,*) 'cfl      = ', cfl
   write(*,*) 'no cells = ', ncel
   write(*,*) 'flux     = ', trim(flux)
   write(*,*) 'recon    = ', trim(recon)
   write(*,*) 'lim var  = ', trim(var)
   write(*,*) 'efix     = ', efix

end subroutine read_input

! Main function where execution starts
program main
   use constants
   use TestData
   implicit none
   real,allocatable  :: u0(:,:), u(:,:), res(:,:), xc(:)
   real              :: dx, dt, t
   integer           :: i, iter
   real,external     :: compute_dt

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
