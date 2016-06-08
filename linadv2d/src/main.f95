program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: co0, co1, res

   integer :: fid

   call read_input

   ! file id for saving solution
   fileid_sol = 0

   ng     = 3
   nrk    = 3

   if(nrk.eq.2)then
      ark(1) = 0.0
      ark(2) = 0.5
   else if(nrk.eq.3)then
      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0
   endif

   allocate( co0((nx+2*ng)*(ny+2*ng)) )
   allocate( co1((nx+2*ng)*(ny+2*ng)) )
   allocate( res((nx+2)*(ny+2)) )

   call solveFVM(co0, co1, res)

end program main
