module comvar
   implicit none

   integer :: nx, ny, ng
   real    :: xmin, xmax, ymin, ymax, dx, dy, dt
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time

   integer :: nrk
   real    :: ark(3)

   real,parameter :: M_PI = 4.0*atan(1.0)

   integer :: fileid_sol

   integer :: fluxtype
   integer :: iupwind=1

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: no=0, yes=1

   integer :: xperiod, yperiod

   contains

      subroutine wave_speed(x, y, speed)
      implicit none
      real :: x, y, speed(2)

      speed(1) = -y
      speed(2) =  x
      end subroutine wave_speed

end module comvar
