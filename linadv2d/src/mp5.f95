real function minmod2(a, b)
   implicit none

   real :: a, b

   if(a*b > 0.0)then
      minmod2 = sign(1.0, a) * min(abs(a), abs(b))
   else
      minmod2 = 0.0
   endif

end function minmod2

real function minmod3(a, b, c)
   implicit none

   real :: a, b, c

   if(a*b > 0.0 .and. b*c > 0.0)then
      minmod3 = sign(1.0, a) * min(abs(a), abs(b), abs(c))
   else
      minmod3 = 0.0
   endif

end function minmod3

real function minmod4(a, b, c, d)
   implicit none

   real :: a, b, c, d

   if(a*b > 0.0 .and. b*c > 0.0 .and. c*d > 0.0)then
      minmod4 = sign(1.0, a) * min(abs(a), abs(b), abs(c), abs(d))
   else
      minmod4 = 0.0
   endif

end function minmod4

real function median(a,b,c)
    implicit none
    real :: a, b, c

    real :: minmod2

    median = a + minmod2(b-a, c-a)
end function median

! Scheme of Suresh and Huynh
subroutine mp5(um2,um1,u0,up1,up2,u)
   implicit none
   real :: um2, um1, u0, up1, up2, u
   real, parameter :: eps = 1.0e-13, alpha = 4.0
   real :: ump, d0, dm1, dp1, dlm4, drm4, uul, uav, umd, ulc
   real :: umin, umax, minmod2, minmod3, minmod4, median

   u = (2.0*um2 - 13.0*um1 + 47.0*u0 + 27.0*up1 - 3.0*up2)/60.0
   ump = u0 + minmod2(up1-u0, alpha*(u0-um1))
   if ((u - u0)*(u - ump) < eps) then
      return
   endif
   d0 = um1 + up1 - 2.0*u0
   dm1= um2 + u0  - 2.0*um1
   dp1= u0  + up2 - 2.0*up1
   dlm4 = minmod4(4*dm1 - d0, 4*d0-dm1, dm1, d0)
   drm4 = minmod4(4*d0 - dp1, 4*dp1-d0, d0,  dp1)
   uul = u0 + alpha*(u0 - um1)
   uav = 0.5*(u0 + up1)
   umd = uav - 0.5*drm4
   ulc = u0 + 0.5*(u0 - um1) + (4.0/3.0)*dlm4
   umin = max(min(u0,up1,umd), min(u0,uul,ulc))
   umax = min(max(u0,up1,umd), max(u0,uul,ulc))
   u = median(u, umin, umax)
end subroutine mp5
