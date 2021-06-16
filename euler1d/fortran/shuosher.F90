module TestData
   use constants
   implicit none
   real,parameter :: xmin = -5.0, xmax = 5.0
   real,parameter :: gam = 1.4
   real,parameter :: Tf = 1.8
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

   if(x < -4.0)then
      v(1) = 3.857143
      v(2) = 2.699369
      v(3) = 10.33333
   else
      v(1) = 1.0 + 0.2*sin(5.0*x)
      v(2) = 0.0
      v(3) = 1.0
   endif

   call prim2con(v, u)
end subroutine initial_condition
