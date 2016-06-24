      subroutine LVQ_method_4(co,Fx,Gy)
      use comvar

      implicit none
      integer:: i, j, ii, jj
      real   :: R, S, lambda
      real, intent (inout)  :: Fx(-ng+1:nx+ng, -ng+1:ny+ng),Gy(-ng+1:nx+ng, -ng+1:ny+ng)
      real, intent (in)  :: co(-ng+1:nx+ng, -ng+1:ny+ng)
      real :: u, v, xf, yf, speed(2) 

     call LVQ_method_3(co,Fx,Gy)
 
      lambda = dt/dx
      do i=1,nx+1
        do j=1,ny
           xf = xmin + (i-1)*dx
           yf = ymin + (j-1)*dy + 0.5*dy
           call wave_speed(xf, yf, speed)
           u = speed(1)
           v = speed(2)

           if(v .gt. 0.0)then 
              jj = j+1
           else
              jj = j
           endif
           R = co(i,j)-co(i-1,j)
           S = 0.5*abs(u)*(1.0-lambda*abs(u))*R
           Gy(i,jj) = Gy(i,jj) + lambda*v*S
           Gy(i-1,jj) = Gy(i-1,jj) - lambda*v*S
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

           if(u .gt. 0.0)then
              ii = i+1
           else
              ii = i
           endif
           R = co(i,j+1)-co(i,j)
           S = 0.5*abs(v)*(1.0-lambda*abs(v))*R
           Fx(ii,j) = Fx(ii,j) + lambda*u*S
           Fx(ii,j-1) = Fx(ii,j-1) - lambda*u*S
        enddo
      enddo

      end subroutine LVQ_method_4
