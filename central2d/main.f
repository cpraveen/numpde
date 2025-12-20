      program main
      include 'dims.h'
      real u(nxd,nyd,mn)

      AMATH_PI = 4.0*atan(1.0)

c     Size of computational domain
      xmin  = -5.0
      xmax  = +5.0
      ymin  = -5.0
      ymax  = +5.0

c     Grid spacing
      dx    = (xmax - xmin)/nx
      dy    = (ymax - ymin)/ny

c     Ratio of specific heats
      gamma = 1.4

c     CFL number
      cfl   = 0.45

c     Limiter
      theta = 2

c     Final time
      tf    = 500.0

c     Vortex parameters
      x0   = 0.0
      y0   = 0.0
      beta = 5.0

      print*,'nx, ny =', nx, ny
      print*,'dx, dy =', dx, dy
      print*,'xmin, xmax =', xmin, xmax
      print*,'ymin, ymax =', ymin, ymax
      print*,'cfl =',cfl
      print*,'tf  =', tf
      print*,'theta =', theta

c     Set initial condition
      print*,'Setting initial conditions'
      do i = md+1, nx+md
         do j = md+1, ny+md
            x    = xmin + (i-md-1)*dx + 0.5*dx
            y    = ymin + (j-md-1)*dy + 0.5*dy
            r    = sqrt( (x-x0)**2 + (y-y0)**2 )
            f1   = (0.5*beta/AMATH_PI)*exp( 0.5*(1.0 - r*r) )
            dvx  = -(y-y0)*f1
            dvy  = +(x-x0)*f1
            drho = 1.0 -
     1           (gamma-1.0)*beta*beta*exp(1.0-r*r)/
     2           (8.0*gamma*AMATH_PI**2)
            drho = drho**(1.0/(gamma-1.0))
            dp   = drho**gamma

            rho = drho
            vx  = dvx
            vy  = dvy
            p   = dp

            u(i,j,1) = rho
            u(i,j,2) = rho*vx
            u(i,j,3) = rho*vy
            u(i,j,4) = p/(gamma-1.0) + 0.5*rho*(vx**2 + vy**2)
         enddo
      enddo

      call savesol(gamma, xmin, ymin, dx, dy, u)
      call system("cp -f out.vtk 0.vtk")
      call system("cp -f line.dat line0.dat")
c     stop "Early stop"

c     Solve until time = tf
      print*,'Starting flow solution'
      call EULER2D(dx,dy,cfl,gamma,theta,tf,u)

      call savesol(gamma, xmin, ymin, dx, dy, u)

      stop
      end
