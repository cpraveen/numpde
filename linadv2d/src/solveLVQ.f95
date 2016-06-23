subroutine solveLVQ(co0,co1,res)

   use comvar

   implicit none

   real :: co0(-ng+1:nx+ng, -ng+1:ny+ng)
   real :: co1(-ng+1:nx+ng, -ng+1:ny+ng)
   real :: res(0:nx+1, 0:ny+1)

   integer :: it, i, j, rks
   real    :: xf, yf
   real    :: xflux, yflux, speed(2)
   real    :: time, lambda
   logical :: tostop

   integer :: ii, jj
   real    :: R1, R2, S, r
   real    :: F(-ng+1:nx+ng, -ng+1:ny+ng), G(-ng+1:nx+ng, -ng+1:ny+ng)

   ! set initial condition
   call init_cond(co1)
   call periodic(co1)
   call savesol(0.0, co1)
   call timestep()

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif
      lambda = dt/dx
      co0 = co1
      res = 0.0; F = 0.0; G = 0.0
 
      if(mthd == method1) call LVQ_method_1(co0,F,G)
      if(mthd == method2) call LVQ_method_2(co0,F,G)
      if(mthd == method3) call LVQ_method_3(co0,F,G)
      if(mthd == method4) call LVQ_method_4(co0,F,G)
      if(mthd == method5) call LVQ_method_5(co0,F,G)
      if(mthd == method6) call LVQ_method_6(co0,F,G)

      ! update conserved variables
      do i=1,nx
        do j=1,ny
           res(i,j) = (F(i+1,j) - F(i,j) + G(i,j+1) - G(i,j))
           co1(i,j) = co0(i,j) - lambda*res(i,j)
        enddo
      enddo
      call periodic(co1)
      it = it + 1
      time = time + dt
      write(*,'(I6,F10.2)')it,time

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call savesol(time, co1)
      endif

   enddo ! time iteration loop

   call error(time, co1)
end subroutine solveLVQ
