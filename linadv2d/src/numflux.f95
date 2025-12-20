subroutine numflux(lx, ly, speed, conl, conr, flux)
   use comvar
   implicit none

   real :: lx, ly, speed(2), conl, conr, flux

   real :: vn

   vn = lx*speed(1) + ly*speed(2)
   if(vn.gt.0.0)then
      flux = vn * conl
   else
      flux = vn * conr
   endif

end subroutine numflux
