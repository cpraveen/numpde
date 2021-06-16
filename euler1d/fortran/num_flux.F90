! Compute Lax-Friedrich flux
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

! Numerical flux of Roe
! ul, ur : conserved variables
subroutine roe_flux(ul, ur, nflux)
   use constants
   use TestData, only : gam
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), fl(nvar), fr(nvar), lam(nvar), HL, HR
   real :: alpha(nvar), r1(nvar), r2(nvar), r3(nvar)
   real ::  u, H, a, RT, du(nvar)

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   HL = (ul(3) + vl(3))/vl(1)
   HR = (ur(3) + vr(3))/vr(1)

   ! Roe averages
   RT  = sqrt(vr(1)/vl(1))
   u   = (vl(2) + RT*vr(2))/(1.0+RT)
   H   = (HL + RT*HR)/(1.0+RT)
   a   = sqrt((gam-1.0)*(H - 0.5*u**2))

   lam(1) = u - a
   lam(2) = u
   lam(3) = u + a

   du = ur - ul

   alpha(2) = (gam - 1.0) * ((H-u**2)*du(1) + u*du(2) - du(3)) / a**2
   alpha(1) = 0.5 * ((u+a)*du(1) - du(2) - a*alpha(2)) / a
   alpha(3) = du(1) - alpha(1) - alpha(2)

   r1(1) = 1.0
   r1(2) = u - a
   r1(3) = H - u*a

   r2(1) = 1.0
   r2(2) = u
   r2(3) = 0.5*u**2

   r3(1) = 1.0
   r3(2) = u + a
   r3(3) = H + u*a

   call euler_flux(vl, fl)
   call euler_flux(vr, fr)

   nflux = 0.5*(fl + fr) &
           - 0.5*(alpha(1)*abs(lam(1))*r1 + alpha(2)*abs(lam(2))*r2 + alpha(3)*abs(lam(3))*r3)
end subroutine roe_flux

! Minimum and maximum wave speed estimates for Riemann data
! vl, vr : primitive variables
subroutine hll_speed(vl, vr, sl, sr)
   use constants
   use TestData, only : gam
   implicit none
   real,intent(in)    :: vl(nvar), vr(nvar)
   real,intent(inout) :: sl, sr
   ! Local variables
   real :: HL, HR, laml, lamr
   real :: u, H, a, upa, uma, RT

   HL   = gam*vl(3)/((gam-1.0)*vl(1)) + 0.5*vl(2)**2
   laml = abs(vl(2)) - sqrt(gam*vl(3)/vl(1))

   HR   = gam*vr(3)/((gam-1.0)*vr(1)) + 0.5*vr(2)**2
   lamr = abs(vr(2)) + sqrt(gam*vr(3)/vr(1))

   ! Roe average speed
   RT  = sqrt(vr(1)/vl(1))
   u   = (vl(2) + RT*vr(2))/(1.0+RT)
   H   = (HL + RT*HR)/(1.0+RT)
   a   = sqrt((gam-1.0)*(H-0.5*u**2))

   uma = u - a
   upa = u + a
   sl  = min(laml, uma)
   sr  = max(lamr, upa)

end subroutine hll_speed

! Numerical flux from HLL solver
! ul, ur : conserved variables
subroutine hll_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), sl, sr, fl(nvar), fr(nvar)

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   call hll_speed(vl, vr, sl, sr)

   if(sl > 0.0) then
      call euler_flux(vl, nflux)
   else if (sr < 0.0) then
      call euler_flux(vr, nflux)
   else
      call euler_flux(vl, fl)
      call euler_flux(vr, fr)
      nflux = (1.0/(sr-sl)) * (sr*fl - sl*fr + sl*sr*(ur - ul))
   endif
end subroutine hll_flux

! Numerical flux from HLLC solver
! ul, ur : conserved variables
subroutine hllc_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), sl, sr, fl(nvar), fr(nvar), ps
   real :: us(nvar), sm

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   call hll_speed(vl, vr, sl, sr)

   if(sl > 0.0)then ! supersonic to right
      call euler_flux(vl, nflux)
   else if(sr < 0.0)then ! supersonic to left
      call euler_flux(vr, nflux)
   else ! subsonic case
      ! u*
      sm =  (ur(2)*(sr-vr(2)) - ul(2)*(sl - vl(2))-(vr(3)-vl(3))) &
          / (vr(1)*(sr -vr(2))- vl(1)*(sl - vl(2)))
      ! p*
      ps = vl(3) + vl(1)*(sl - vl(2))*(sm - vl(2))
      if(0.0 <= sm)then
         call euler_flux(vl, fl)
         us(1) = ul(1) * (sl - vl(2))/(sl - sm)
         us(2) = us(1)*sm
         us(3) = ((sl-vl(2))*ul(3) + ps*sm - vl(3)*vl(2)) / (sl - sm)
         nflux = fl + sl*(us - ul)
      else
         call euler_flux(vr, fr)
         us(1) = ur(1) * (sr - vr(2))/(sr - sm)
         us(2) = us(1)*sm
         us(3) = ((sr-vr(2))*ur(3) + ps*sm - vr(3)*vr(2)) / (sr - sm)
         nflux = fr + sr*(us - ur)
      endif
   endif
end subroutine hllc_flux
