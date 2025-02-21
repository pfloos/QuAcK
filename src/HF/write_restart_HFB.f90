
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
  double precision                   :: val_print

! Dump results

  open(unit=iunit,form='unformatted',file='hfb_bin')
  write(iunit) nBas,nOrb 
  write(iunit) chem_pot
  do iorb=1,nOrb 
   do ibas=1,nBas
    val_print = c(ibas,nOrb+1-iorb)
    if(abs(val_print)<1d-8) val_print=0d0
    write(iunit) val_print
   enddo
  enddo
  do iorb=1,nOrb 
   val_print = Occ(nOrb+1-iorb)
   if(abs(val_print)<1d-8) val_print=0d0
   write(iunit) val_print
  enddo
  write(iunit) iunit 
  close(iunit)

end subroutine 
