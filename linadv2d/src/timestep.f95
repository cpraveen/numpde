subroutine timestep()
   use comvar
   implicit none

   integer :: i, j
   real    :: speed, eig

   speed = 0.0

   do i=1,nx
      do j=1,ny
         speed = 1.0
      enddo
   enddo

   dt = dx/speed
   dt = cfl*dt

   write(*,*)'Time step =', dt

end subroutine timestep
