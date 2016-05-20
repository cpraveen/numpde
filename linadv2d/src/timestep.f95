subroutine timestep()
   use comvar
   implicit none

   integer :: i, j
   real    :: x, y, speed(2)

   dt = 1.0e20
   do i=1,nx
      do j=1,ny
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         call wave_speed(x, y, speed)
         dt = min(dt, 1.0/(abs(speed(1))/dx + abs(speed(2))/dy + 1.0e-20))
      enddo
   enddo

   dt = cfl*dt

   write(*,*)'Time step =', dt

end subroutine timestep
