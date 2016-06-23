      subroutine LVQ_method_6(co,Fx,Gy)
      use comvar

      implicit none
      integer:: i, j, ii, jj
      real   :: R, S, R2, lambda, lambdasq
      real, intent (inout)  :: Fx(-ng+1:nx+ng, -ng+1:ny+ng),Gy(-ng+1:nx+ng, -ng+1:ny+ng)
      real, intent (in)  :: co(-ng+1:nx+ng, -ng+1:ny+ng)
      real :: u, v, xf, yf, speed(2)
      
      call LVQ_method_5(co,Fx,Gy)

      lambda = dt/dx; lambdasq = lambda*lambda
      do i=1,nx+1
        do j=1,ny
           xf = xmin + (i-1)*dx
           yf = ymin + (j-1)*dy + 0.5*dy
           call wave_speed(xf, yf, speed)
           u = speed(1)
           v = speed(2)

           if(u .lt. 0.0)then 
              ii = i
           else
              ii = i-1 
           endif
           R2 = co(ii+1,j) - 2.0*co(ii,j) + co(ii-1,j)
           S = (1.0/6.0)*u*(1.0-u*u*lambdasq)*R2
           Gy(i-1,j+1) = Gy(i-1,j+1) - lambda*v*S
           if(v .lt. 0.0)then 
              jj = j
           else
              jj = j-1 
           endif
           R2 = co(i,jj+1) - 2.0*co(i,jj) + co(i,jj-1)
           S = (1.0/6.0)*u*(1.0-lambdasq*u*u)*R2
           Gy(i,j+1) = Gy(i,j+1) + lambda*v*S
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
           if(v .lt. 0.0)then 
              jj = j
           else
              jj = j-1 
           endif
           R2 = co(i,jj+1) - 2.0*co(i,jj) + co(i,jj-1)
           S = (1.0/6.0)*v*(1.0-lambdasq*v*v)*R2
           Fx(i+1,j-1) = Fx(i+1,j-1) - lambda*u*S
           if(u .lt. 0.0)then 
              ii = i
           else
              ii = i-1 
           endif
           R2 = co(ii+1,j) - 2.0*co(ii,j) + co(ii-1,j)
           S = (1.0/6.0)*v*(1.0-lambdasq*v*v)*R2
           Fx(i+1,j) = Fx(i+1,j) + lambda*u*S
        enddo
      enddo

      end subroutine LVQ_method_6
