
! ---

subroutine write_restart_HFB(nBas, nOrb, Occ, c, chem_pot)

! Write the binary file used to restart calcs

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nOrb
  double precision,intent(in)        :: chem_pot
  double precision,intent(in)        :: Occ(nOrb)
  double precision,intent(in)        :: c(nBas,nOrb)

! Local variables

  integer                            :: ibas,iorb,iunit=666

! Dump results

  open(unit=iunit,form='unformatted',file='hfb_bin')
  write(iunit) nBas,nOrb 
  write(iunit) chem_pot
  do iorb=1,nOrb 
   write(iunit) c(:,nOrb+1-iorb) 
  enddo
  do iorb=1,nOrb 
   write(iunit) Occ(nOrb+1-iorb) 
  enddo
  write(iunit) iunit 
  close(iunit)

end subroutine 
