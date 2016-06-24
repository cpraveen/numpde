subroutine read_input
   use comvar
   implicit none
   integer :: fid
   character(40) :: fluxtype_s, limtype_s, mthd_s

   open(10,file='inp.dat')
   read(10,*) nx
   read(10,*) ny
   read(10,*) itmax
   read(10,*) itsave
   read(10,*) fluxtype_s
   read(10,*) limtype_s
   read(10,*) mthd_s
   read(10,*) cfl
   read(10,*) testcase

   print*,'nx, ny =', nx, ny
   print*,'cfl    =', cfl
   print*,'Flux type is ', fluxtype_s
   if(fluxtype_s == "upwind")then
      fluxtype = iupwind
   else if(fluxtype_s == "mda")then
      fluxtype = imda
   else
      print*,'Unknown fluxtype =', fluxtype
      stop
   endif

   print*,'Reconstruction type is ', limtype_s
   if(limtype_s == "first")then
      limtype = iford
   else if(limtype_s == "muscl3")then
      limtype = imuscl3
   else if(limtype_s == "minmod")then
      limtype = immod
   else if(limtype_s == "weno5")then
      limtype = iweno5
   else if(limtype_s == "mp5")then
      limtype = imp5
   else
      print*,'Unknown limtype =', limtype_s
      stop
   endif

   print*,'LeVeque Method type is ', mthd_s
   if(mthd_s == "method1")then
      mthd = method1
   else if(mthd_s == "method2")then
      mthd = method2
   else if(mthd_s == "method3")then
      mthd = method3
   else if(mthd_s == "method4")then
      mthd = method4
   else if(mthd_s == "method5")then
      mthd = method5
   else if(mthd_s == "method6")then
      mthd = method6
   else
      print*,'Unknown method =', mthd_s
      stop
   endif
end subroutine read_input
