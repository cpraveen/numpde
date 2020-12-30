subroutine timestep()
   use comvar
   implicit none
   ! Local variables
   integer :: i, j
   real    :: x, y, speed(2), dt1, dt2

   ! If dt > 0 is given in input file, use it
   if(dt .gt. 0.0)then
      write(*,*)'Using time step from input file'
      write(*,*)'Time step =', dt
      return
   endif

   dt1 = 1.0e20
   dt2 = 1.0e20
   do i=1,nx
      do j=1,ny
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         call wave_speed(x, y, speed)
         dt1 = min(dt1, 1.0/(abs(speed(1))/dx + abs(speed(2))/dy + 1.0e-20))
         dt2 = min(dt2, dx/(abs(speed(1)) + 1.0e-20))
         dt2 = min(dt2, dy/(abs(speed(2)) + 1.0e-20))
      enddo
   enddo

   if(fluxtype == iupwind)then
      dt = dt1
   else if(fluxtype == imda)then
      dt = dt2
   else if(fluxtype == ilw)then
      dt = 0.72 * dt1
   else
      print*,'timestep: unknown fluxtype =',fluxtype
      stop
   endif

   dt = cfl*dt

   write(*,*)'Time step =', dt

end subroutine timestep
