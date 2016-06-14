subroutine reconstruct(conjm2, conjm1, conj, conjp1, conjp2, conl)
   use comvar
   implicit none

   real    :: conjm2, conjm1, conj, conjp1, conjp2
   real    :: conl

   integer :: i
   real    :: minmod
   real,parameter :: kkk=1.0/3.0
   real,parameter :: beta = 2.0

   ! reconstructed states
   if(limtype == iford)then
   ! first order
      conl = conj
   else if(limtype == imuscl3)then
   !muscl scheme
      conl = conj   + 0.25*( (1.0-kkk)*(conj - conjm1) &
                                 + (1.0+kkk)*(conjp1 - conj) )
   else if(limtype == immod)then
   ! minmod limiter
      conl = conj + 0.5*minmod( beta*(conj-conjm1), &
                                0.5*(conjp1-conjm1), &
                                beta*(conjp1-conj) )
   else if(limtype == iweno5)then
      call weno5(conjm2, conjm1, conj, conjp1, conjp2, conl)
   else if(limtype == imp5)then
      call mp5(conjm2, conjm1, conj, conjp1, conjp2, conl)
   endif

end subroutine reconstruct
