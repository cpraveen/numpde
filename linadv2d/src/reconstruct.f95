subroutine reconstruct(conjm2, conjm1, conj, conjp1, conjp2, conl)
   use comvar
   implicit none

   real    :: conjm2, conjm1, conj, conjp1, conjp2
   real    :: conl

   integer :: i
   real    :: kkk=1.0/3.0
   real    :: minmod

   ! reconstructed states
   if(limtype == ford)then
   ! first order
      conl = conj
   else if(limtype == muscl3)then
   !muscl scheme
      conl = conj   + 0.25*( (1.0-kkk)*(conj - conjm1) &
                                 + (1.0+kkk)*(conjp1 - conj) )
   else if(limtype == mmod)then
   ! minmod limiter
      conl = conj + 0.5*minmod( conj-conjm1, &
                                0.5*(conjp1-conjm1), &
                                conjp1-conj )
   endif

end subroutine reconstruct
