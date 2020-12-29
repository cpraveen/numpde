subroutine solveLW(co1, res)
   use comvar
   implicit none
   real :: co1(-ng+1:nx+ng, -ng+1:ny+ng)
   real :: res(0:nx+1,0:ny+1)

   integer :: it, i, j
   real    :: lambda, xf, yf
   real    :: xflux, yflux, speed(2)
   real    :: time
   logical :: tostop

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
      lambda = dt/dx/dy

      res = 0.0
         
      ! x fluxes: between (i,j) and (i+1,j)
      do i=0,nx
         do j=1,ny
            xf = xmin + i*dx
            yf = ymin + (j-1)*dy + 0.5*dy
            call wave_speed(xf, yf, speed)
            call lwflux(speed(1), speed(2), dx, dy, dt, &
                        co1(i,j),  co1(i+1,j), &
                        co1(i,j-1),co1(i+1,j-1), &
                        co1(i,j+1),co1(i+1,j+1), &
                        xflux)
            res(i,j)   = res(i,j)   + dy*xflux
            res(i+1,j) = res(i+1,j) - dy*xflux
         enddo
      enddo
         
      ! y fluxes: between (i,j) and (i,j+1)
      do j=0,ny
         do i=1,nx
            xf = xmin + (i-1)*dx + 0.5*dx
            yf = ymin + j*dy
            call wave_speed(xf, yf, speed)
            call lwflux(speed(2), -speed(1), dy, dx, dt, &
                        co1(i,j),  co1(i,j+1), &
                        co1(i+1,j),co1(i+1,j+1), &
                        co1(i-1,j),co1(i-1,j+1), &
                        yflux)
            res(i,j)   = res(i,j)   + dx*yflux
            res(i,j+1) = res(i,j+1) - dx*yflux
         enddo
      enddo

      ! update conserved variables
      do i=1,nx
         do j=1,ny
            co1(i,j) = co1(i,j) - lambda*res(i,j)
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

end subroutine solveLW

! LW flux in x-direction for face between (ql,qr)
!          -----------
!           qul | qur
!          -----------
!           ql  | qr
!          -----------
!           qdl | qdr
!          -----------
subroutine lwflux(u,v,dx,dy,dt,ql,qr,qdl,qdr,qul,qur,flux)
   implicit none
   real,intent(in)    :: u,v,dx,dy,dt,ql,qr,qdl,qdr,qul,qur
   real,intent(inout) :: flux
   ! Local variables
   real :: qt

   qt = - u * (qr - ql)/dx - v * (qul - qdl + qur - qdr)/(4.0*dy)
   flux = u * (ql + qr)/2.0  + 0.5 * dt * u * qt

end subroutine lwflux
