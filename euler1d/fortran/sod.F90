module TestData
   use constants
   implicit none
   real,parameter :: xmin = 0.0, xmax = 1.0
   real,parameter :: gam = 1.4
   real,parameter :: Tf = 0.2
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

   if(x < 0.5)then
      v(1) = 1.0
      v(2) = 0.0
      v(3) = 1.0
   else
      v(1) = 0.125
      v(2) = 0.0
      v(3) = 0.1
   endif

   call prim2con(v, u)
end subroutine initial_condition
