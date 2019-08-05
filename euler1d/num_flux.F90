subroutine roe_flux(ul, ur, nflux)
   use constants
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

subroutine hll_speed(ul, ur, sl, sr)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: sl, sr
   ! Local variables
   real :: vl(nvar), vr(nvar), HL, HR, laml, lamr 
   real ::  u, H, a, upa, uma, RT

   call con2prim(ul, vl)
   call con2prim(ur, vr)
   
   HL   = (ul(3) + vl(3))/vl(1)
   laml = abs(vl(2)) - sqrt(gam*vl(3)/vl(1))
   
   HR   = (ur(3) + vr(3))/vr(1)
   lamr = abs(vr(2)) + sqrt(gam*vr(3)/vr(1))

   RT  = sqrt(vr(1)/vl(1))
   u   = (vl(2) + RT*vr(2))/(1.0+RT)
   H   = (HL + RT*HR)/(1.0+RT)
   a   = sqrt((gam-1.0)*(H-0.5*u**2))

   uma = u - a
   upa = u + a
   sl  = min(laml, uma)
   sr  = max(lamr, upa)
   
end subroutine hll_speed

subroutine hll_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), sl, sr, fl(nvar), fr(nvar)

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   call hll_speed(ul, ur, sl, sr)
   
   if(sl > 0.0) then
      call euler_flux(vl, nflux)
   else if (sr < 0.0) then
      call euler_flux(vr, nflux)
   else 
      call euler_flux(vl, fl)
      call euler_flux(vr, fr)
      nflux = (1.0/(sr-sl)) * ( sr*fl - sl*fr + sl*sr*(ur - ul))
   endif
end subroutine hll_flux

subroutine hllc_flux(ul, ur, nflux)
   use constants
   implicit none
   real,intent(in)    :: ul(nvar), ur(nvar)
   real,intent(inout) :: nflux(nvar)
   ! Local variables
   real :: vl(nvar), vr(nvar), sl, sr, fl(nvar), fr(nvar)
   real :: us(nvar), vs(nvar), fa(nvar), sm, sa, ua(nvar), va(nvar), usa(nvar)

   call con2prim(ul, vl)
   call con2prim(ur, vr)

   call hll_speed(ul, ur, sl, sr)
   
   if(sl > 0.0)then ! supersonic to right
      call euler_flux(vl, nflux)
   else if(sr < 0.0)then ! supersonic to left
      call euler_flux(vr, nflux)
   else ! subsonic case
      call euler_flux(vl, fl)
      call euler_flux(vr, fr)
      ! HLL state
      us = (sr*ur - sl*ul - (fr - fl))/(sr - sl)
      call con2prim(us, vs)
      sm = vs(2)
      if(0.0 <= sm)then
         fa = fl
         ua = ul
         sa = sl
      else
         fa = fr
         ua = ur
         sa = sr
      endif

      call con2prim(ua, va)
      usa(1) = (sa - va(2))/(sa - sm) * ua(1)
      usa(2) = usa(1)*sm
      usa(3) = ((sa-va(2))*ua(3) + vs(3)*sm - va(3)*va(2)) / (sa - sm)

      nflux = fa + sa*(usa - ua)
   endif
end subroutine hllc_flux