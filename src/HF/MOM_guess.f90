subroutine MOM_guess(nO, nBas, nOrb, occ,c,cGuess,eHF)

! Prepares the MOM guess based on the gorund state orbitals c and occupations, i.e.
! It resorts the columns of c such that the occupied orbitals are the first nO

  implicit none
  include 'parameters.h'

! Input variables
  integer, intent(in)             :: nO,nBas, nOrb
  integer, intent(in)             :: occ(nO)
  double precision, intent(in)    :: c(nBas,nOrb)

! Input/Output variables
  double precision, intent(out)   :: cGuess(nBas, nOrb)
  double precision, intent(out)   :: eHF(nOrb)

! Local variables
  integer                         :: i,p,k
  double precision,allocatable    :: etmp(:)

  allocate(etmp(nOrb))

! Occupied orbitals
  do i = 1, nO
     if(occ(i) == 0) then
       print *, "Number of occupied orbitals is too small !"
       stop
     end if
     cGuess(:,i) = c(:,occ(i))
     etmp(i) = eHF(occ(i))
  end do

! Virtual orbitals
  k = 0
  do p = 1, nOrb
    if(.not. any(occ == p)) then 
      k = k + 1
      cGuess(:,nO + k) = c(:,p)
      etmp(nO+k) = eHF(p)
    end if
  end do

  eHF(:) = etmp(:)
  
  deallocate(etmp)
  
  if(nO + k  /= nOrb) then
    print *, "Number of occupied orbitals is not correct !"
    stop
  end if

end subroutine
