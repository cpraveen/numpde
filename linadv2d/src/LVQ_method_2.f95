      subroutine LVQ_method_2(co,Fx,Gy)
      use comvar

      implicit none
      integer:: i,j,ii,jj
      real   :: R, lambda
      real, intent (inout)  :: Fx(-ng+1:nx+ng, -ng+1:ny+ng)
      real, intent (inout)  :: Gy(-ng+1:nx+ng, -ng+1:ny+ng)
      real, intent (in)  :: co(-ng+1:nx+ng, -ng+1:ny+ng)
      real :: u, v, xf, yf, speed(2)

      call LVQ_method_1(co,Fx,Gy)

      lambda = dt/dx
      do i=1,nx+1
        do j=1,ny
           xf = xmin + (i-1)*dx
           yf = ymin + (j-1)*dy + 0.5*dy
           call wave_speed(xf, yf, speed)
           u = speed(1)
           v = speed(2)

           R = co(i,j) - co(i-1,j)
           if(u .gt. 0.0)then 
             ii = i
           else
             ii = i-1
           endif
           if(v .gt. 0.0)then 
             jj = j+1
           else
             jj = j
           endif
           Gy(ii,jj)   = Gy(ii,jj) - 0.5*lambda*u*v*R
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

          R = co(i,j) - co(i,j-1)

          if(v .gt. 0.0)then 
             jj = j
          else
             jj = j-1
          endif
          if(u .gt. 0.0)then 
             ii = i+1
          else
             ii = i
          endif
          Fx(ii,jj)   = Fx(ii,jj) - 0.5*lambda*u*v*R
        enddo
      enddo

      end subroutine LVQ_method_2
