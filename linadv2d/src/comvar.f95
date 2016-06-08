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
   integer :: iupwind=1,imda=2

   integer :: limtype
   integer,parameter :: iford=0, imuscl3=1, immod=2, iweno5=3

   integer :: testcase

   integer,parameter :: no=0, yes=1

   integer :: xperiod, yperiod

   contains

!     wave speed in the advection equation
      subroutine wave_speed(x, y, speed)
      implicit none
      real :: x, y, speed(2)

      speed(1) = -y
      speed(2) =  x
      end subroutine wave_speed

end module comvar
