subroutine non_occupied(nO,nOrb,occupations, virtual)

! Given occupations, this subroutine determines the virtual (not-occupied) orbitals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nO, nOrb
  integer,intent(in)            :: occupations(nO)

! Local variables
  integer                       :: p, counterVirtual

! Output variables

  integer,intent(out)           :: virtual(nOrb - nO)

  counterVirtual  = 0
  do p = 1,nOrb
    if( .not. any(p == occupations(:))) then
      
      counterVirtual = counterVirtual + 1
      virtual(counterVirtual) = p
    
    end if
  end do 
end subroutine
