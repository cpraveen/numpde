subroutine getfilename(filename, fileid)
   implicit none

   character(len=*) :: filename
   integer :: fileid

   if(fileid <= 9)then
      write(unit=filename, fmt='(A5,I1)') trim(filename)//'00', fileid
   else if(fileid <= 99)then
      write(unit=filename, fmt='(A4,I2)') trim(filename)//'0', fileid
   else if(fileid <= 999)then
      write(unit=filename, fmt='(A3,I3)') trim(filename), fileid
   endif

   filename = trim(filename)//'.plt'

   fileid = fileid + 1

end subroutine getfilename
