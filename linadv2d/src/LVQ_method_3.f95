      subroutine LVQ_method_3(co,Fx,Gy)
      use comvar

      implicit none
      integer:: i, j
      real   :: R, S
      real, intent (inout)  :: Fx(-ng+1:nx+ng, -ng+1:ny+ng),Gy(-ng+1:nx+ng, -ng+1:ny+ng)
      real, intent (in)  :: co(-ng+1:nx+ng, -ng+1:ny+ng)
      real :: u, v, xf, yf, speed(2) 

      call LVQ_method_2(co,Fx,Gy)

      ! x fluxes
      do i=1,nx+1
        do j=1,ny
           xf = xmin + (i-1)*dx
           yf = ymin + (j-1)*dy + 0.5*dy
           call wave_speed(xf, yf, speed)
           u = speed(1)
           v = speed(2)

           R = co(i,j)-co(i-1,j)
           S = 0.5*abs(u)*(1.0-(dt/dx)*abs(u))*R
           Fx(i,j) = Fx(i,j) + S
        enddo
      enddo

      ! y fluxes
      do j=1,ny+1
        do i=1,nx
           xf = xmin + (i-1)*dx + 0.5*dx
           yf = ymin + (j-1)*dy
           call wave_speed(xf, yf, speed)
           u = speed(1)
           v = speed(2)

           R = co(i,j)-co(i,j-1)
           S = 0.5*abs(v)*(1.0-(dt/dx)*(abs(v)))*R 
           Gy(i,j) = Gy(i,j) + S
        enddo
      enddo
      end subroutine LVQ_method_3
