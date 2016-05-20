!---------------------------------------------------------------------
! Weno5 reconstruction
!---------------------------------------------------------------------
      subroutine weno5(um2,um1,u0,up1,up2,u)
      implicit none
      real :: um2, um1, u0, up1, up2, u
      real,parameter :: eps = 1.0e-6, gamma1=1.0/10.0, gamma2=3.0/5.0, &
                        gamma3=3.0/10.0
      real :: beta1, beta2, beta3
      real :: u1, u2, u3;
      real :: w1, w2, w3;

      beta1 = (13.0/12.0)*(um2 - 2.0*um1 + u0)**2 + &
        (1.0/4.0)*(um2 - 4.0*um1 + 3.0*u0)**2
      beta2 = (13.0/12.0)*(um1 - 2.0*u0 + up1)**2 + &
        (1.0/4.0)*(um1 - up1)**2
      beta3 = (13.0/12.0)*(u0 - 2.0*up1 + up2)**2 + &
        (1.0/4.0)*(3.0*u0 - 4.0*up1 + up2)**2

      w1 = gamma1 / (eps+beta1)**2
      w2 = gamma2 / (eps+beta2)**2
      w3 = gamma3 / (eps+beta3)**2

      u1 = (1.0/3.0)*um2 - (7.0/6.0)*um1 + (11.0/6.0)*u0
      u2 = -(1.0/6.0)*um1 + (5.0/6.0)*u0 + (1.0/3.0)*up1
      u3 = (1.0/3.0)*u0 + (5.0/6.0)*up1 - (1.0/6.0)*up2

      u = (w1 * u1 + w2 * u2 + w3 * u3)/(w1 + w2 + w3)

      return
      end subroutine weno5
