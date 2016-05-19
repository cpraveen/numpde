subroutine savesol(t, co)
   use comvar
   implicit none

   real    :: t
   real    :: co(-ng+1:nx+ng, -ng+1:ny+ng)

   integer :: i, j
   real    :: x, y
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "Linear Advection"'
   write(10,*)'VARIABLES = "x", "y", "Sol"'
   write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=',t,', I=',nx,', J=',ny,&
              ', DATAPACKING=POINT'

   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         write(10,'(3E24.14)') x, y, co(i,j)
      enddo
   enddo

   close(10)

end subroutine savesol
