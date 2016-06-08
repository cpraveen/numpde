subroutine init_cond(co1)
   use comvar
   implicit none

   real    :: co1(-ng+1:nx+ng, -ng+1:ny+ng)

   real    :: x, y
   integer :: i, j

   print*,'Setting initial condition'

   final_time = 2.0*M_PI

   ! periodicity conditions
   xperiod = yes
   yperiod = yes

   xmin = -1.0
   xmax =  1.0
   ymin = -1.0
   ymax =  1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         co1(i,j) = 1.0 + exp(-50.0*((x-0.5)**2 + y**2))
      enddo
   enddo

end subroutine init_cond
