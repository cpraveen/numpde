subroutine error(t, co)
   use comvar
   implicit none

   real    :: t
   real    :: co(-ng+1:nx+ng, -ng+1:ny+ng)

   integer :: i, j
   real    :: du, err_l1, err_l2, err_inf, area
   real    :: co0(-ng+1:nx+ng, -ng+1:ny+ng)

   ! Exact solution if same as initial condition
   call init_cond(co0)

   err_l1  = 0.0
   err_l2  = 0.0
   err_inf = -1.0e20

   do j=1,ny
      do i=1,nx
         du = abs(co(i,j) - co0(i,j))
         err_l1 = err_l1 + du * dx * dy
         err_l2 = err_l2 + du**2 * dx * dy
         err_inf= max(err_inf, du)
      enddo
   enddo
   area   = nx * ny * dx * dy  ! area of domain
   err_l1 = err_l1/area
   err_l2 = sqrt(err_l2/area)

   write(*,*)'Error at time =', t
   write(*,'(2x,"nx",3x,"ny",8x,"dx",12x,"dy",12x,"l1",12x,"l2",11x,"linf")')
   write(*,'(2i5,5e14.6)') nx, ny, dx, dy, err_l1, err_l2, err_inf

end subroutine error
