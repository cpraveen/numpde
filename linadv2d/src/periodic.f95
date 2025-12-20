subroutine periodic(con)
   use comvar
   implicit none

   real :: con(-ng+1:nx+ng,-ng+1:ny+ng)

   integer :: i, j, k, l

   if(xperiod == yes)then
      do j=1,ny
         do k=1,ng
            con(1-k, j) = con(nx-k,j)
            con(nx+k,j) = con(k+1, j)
         enddo   
      enddo
   endif

   if(yperiod == yes)then
      do i=1,nx
         do k=1,ng
            con(i,1-k)  = con(i,ny-k)
            con(i,ny+k) = con(i,k+1)
         enddo   
      enddo
   endif

! We dont use the corner points
    if(xperiod == yes .and. yperiod == yes)then
      do k=1,ng
         do l=1,ng
            con(1-k,1-l)   = con(nx-k,ny-l) ! lower left corner
            con(nx+k,1-l)  = con(k+1,ny-l) ! lower right corner
            con(1-k,ny+l)  = con(nx-k,l+1) ! upper left corner
            con(nx+k,ny+l) = con(k+1,l+1) ! upper right corner
         enddo
      enddo
    endif

end subroutine periodic
