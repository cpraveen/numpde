      subroutine savesol(gamma, xmin, ymin, dx, dy, u)
      include 'dims.h'
      real u(nxd,nyd,mn)

      real om(nxd,nyd)

      character*32 ofile

      om(:,:) = 0.0
      do i=md+1,nx+md
         do j=md+1,ny+md
            v_x = 0.5*( u(i+1,j-1,3) + u(i+1,j+1,3) ) -
     1            0.5*( u(i-1,j-1,3) + u(i-1,j+1,3) )
            v_x = v_x/dx
            u_y = 0.5*( u(i-1,j+1,2) + u(i+1,j+1,2) ) -
     1            0.5*( u(i-1,j-1,2) + u(i+1,j-1,2) )
            u_y = u_y/dy
            om(i,j) = v_x - u_y
         enddo
      enddo

      ofile = 'out.vtk'

      ibeg = md + 1
      iend = nx + md
      jbeg = md + 1
      jend = ny + md

      ni   = iend - ibeg + 1
      nj   = jend - jbeg + 1

      ivtk = 49
      write(*,*)'Writing vtk file ', ofile
      open(unit=ivtk, file=ofile)
      write(ivtk,'("# vtk DataFile Version 2.0")')
      write(ivtk,'("VORTEX")')
      write(ivtk,'("ASCII")')
      write(ivtk,'("DATASET STRUCTURED_GRID")')
      write(ivtk,'("DIMENSIONS ",i10,i10,i10)') ni, nj, 1
      write(ivtk,'("POINTS ",i10," float")') ni*nj
      k = 2
      do j=jbeg,jend
         do i=ibeg,iend
            x   = xmin + (i-md-1)*dx + 0.5*dx
            y   = ymin + (j-md-1)*dy + 0.5*dy
            write(ivtk,*) x, y, 0.0
         enddo
      enddo
      write(ivtk,'("POINT_DATA",i10)') ni*nj
c     Write density
      write(ivtk,'("SCALARS density float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) u(i,j,1)
         enddo
      enddo
c     Write u velocity
      write(ivtk,'("SCALARS vx float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) u(i,j,2)/u(i,j,1)
         enddo
      enddo
c     Write v velocity
      write(ivtk,'("SCALARS vy float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) u(i,j,3)/u(i,j,1)
         enddo
      enddo
c     Write pressure
      write(ivtk,'("SCALARS pressure float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            p = (gamma-1.0)*( u(i,j,4) - 0.5*(u(i,j,2)**2 + u(i,j,3)**2)
     1          / u(i,j,1) )
            write(ivtk,*) p
         enddo
      enddo
c     Write vorticity
      write(ivtk,'("SCALARS vorticity float 1")')
      write(ivtk,'("LOOKUP_TABLE default")')
      do j=jbeg,jend
         do i=ibeg,iend
            write(ivtk,*) om(i,j)
         enddo
      enddo
      close(ivtk)

      ofile = 'line.dat'
      iline = 20
      open(iline, file=ofile)
      j = (ny+2*md)/2
      do i=md+1,nx+md
         x   = xmin + (i-md-1)*dx + 0.5*dx
         y   = ymin + (j-md-1)*dy + 0.5*dy
         p = (gamma-1.0)*( u(i,j,4) - 0.5*(u(i,j,2)**2 + u(i,j,3)**2)
     1        / u(i,j,1) )
         write(iline,*) x, y, u(i,j,2)/u(i,j,1), u(i,j,3)/u(i,j,1), p
      enddo
      close(iline)

      return
      end
