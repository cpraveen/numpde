module TestData
   use constants
   implicit none
   real,parameter :: xmin = 0.0, xmax = 1.0
   real,parameter :: gam = 1.4
   real,parameter :: Tf = 0.012
   integer,parameter :: leftbc = iNeumann, rightbc = iNeumann
end module TestData

! Initial condition as a function of x
subroutine initial_condition(x, u)
   use constants
   implicit none
   real,intent(in)    :: x
   real,intent(inout) :: u(nvar)
   ! Local variables
   real :: v(nvar)

   if(x < 0.8)then
      v(1) = 1.0
      v(2) = -19.59745
      v(3) = 1000.0
   else
      v(1) = 1.0
      v(2) = -19.59745
      v(3) = 0.1
   endif

   call prim2con(v, u)
end subroutine initial_condition
