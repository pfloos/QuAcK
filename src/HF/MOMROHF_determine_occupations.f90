subroutine MOMROHF_determine_occupations(nOa, nOb,nOrb,occupations,                     &
                                         counterSingleA,counterSingleB,counterDouble,   &
                                         counterVirtual,                                &
                                         singlyOccupiedA,singlyOccupiedB,doublyOccupied,&
                                         virtual)

! Determines doubly occupied orbitals and singly occupied orbitals.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOa, nOb, nOrb
  integer,intent(in)            :: occupations(max(nOa,nOb),nspin)

! Local variables
  integer                       :: p

! Output variables

  integer,intent(out)           :: doublyOccupied(max(nOa,nOb))
  integer,intent(out)           :: singlyOccupiedA(nOa),singlyOccupiedB(nOb)
  integer,intent(out)           :: virtual(nOrb - max(nOa,nOb))
  integer,intent(out)           :: counterDouble, counterSingleA, counterSingleB,counterVirtual


! Number of closed, open, and virtual orbitals
  counterDouble = 0
  counterSingleA = 0
  counterSingleB = 0
  
  ! Single alpha and double
  do p=1,nOa
    if(any(occupations(p,1)==occupations(1:nOb,2))) then
      counterDouble  = counterDouble + 1
      doublyOccupied(counterDouble) = occupations(p,1)
    else
      counterSingleA = counterSingleA + 1
      singlyOccupiedA(counterSingleA) = occupations(p,1)
    end if
  end do
  
  ! Single beta
  do p=1,nOb
    if(.not. any(occupations(p,2)==occupations(1:nOa,1))) then
      counterSingleB = counterSingleB + 1
      singlyOccupiedB(counterSingleB) = occupations(p,2)
    end if
  end do

  ! Check occupancy consistence 
  if(     counterSingleA + counterDouble /= nOa &
     .or. counterSingleB + counterDouble /= nOb   ) then
     print *, "Error in occupations : Single and Double occupations do not sum up to number of alpha (resp. beta) electrons"
     stop
  endif
  
  ! Virtuals
  counterVirtual = 0
  do p = 1,nOrb
    if(     (.not. any(p == singlyOccupiedA(1:counterSingleA))) &
      .and. (.not. any(p == singlyOccupiedB(1:counterSingleB))) &
      .and. (.not. any(p ==  doublyOccupied(1:counterDouble)))) then
      
      counterVirtual = counterVirtual + 1 
      virtual(counterVirtual) = p
    
    end if
  end do 
end subroutine
