program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: co0, co1, res

   integer :: fid

   call read_input

   ! file id for saving solution
   fileid_sol = 0

   ng     = 3 ! number of ghost cells
   nrk    = 3 ! number of rk stages, only 2 or 3 is valid

   if(nrk.eq.2)then
      ark(1) = 0.0
      ark(2) = 0.5
   else if(nrk.eq.3)then
      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0
   else
      print*,'Only 2 or 3 stage RK implemented'
      stop
   endif

   allocate( co0((nx+2*ng)*(ny+2*ng)) )
   allocate( co1((nx+2*ng)*(ny+2*ng)) )
   allocate( res((nx+2)*(ny+2)) )

   if(fluxtype==iupwind)then
      call solveFVM(co0, co1, res)
   else if(fluxtype==imda)then
      call solveLVQ(co0, co1, res)
   else
      print*,'Unknown fluxtype'
      stop
   endif
end program main
