
! ---

subroutine write_restart_qsGWB(nBas, nOrb, Occ, c, chem_pot)

! Write the binary file used to restart calcs

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nOrb
  double precision,intent(in)        :: chem_pot
  double precision,intent(inout)     :: Occ(nOrb)
  double precision,intent(inout)     :: c(nBas,nOrb)

! Local variables

  integer                            :: ibas,iorb,iunit=666
  double precision                   :: val_print
  double precision,allocatable       :: Occ_tmp(:)
  double precision,allocatable       :: c_tmp(:,:)

! Dump results

  allocate(Occ_tmp(nOrb),c_tmp(nBas,nOrb))

  open(unit=iunit,form='unformatted',file='qsgwb_bin')
  write(iunit) nBas,nOrb 
  write(iunit) chem_pot
  do iorb=1,nOrb 
   c_tmp(:,iorb) = c(:,nOrb+1-iorb)
   do ibas=1,nBas
    val_print = c(ibas,nOrb+1-iorb)
    if(abs(val_print)<1d-8) val_print=0d0
    write(iunit) val_print
   enddo
  enddo
  do iorb=1,nOrb
   Occ_tmp(iorb) = Occ(nOrb+1-iorb) 
   val_print = Occ(nOrb+1-iorb)
   if(abs(val_print)<1d-8) val_print=0d0
   write(iunit) val_print
  enddo
  write(iunit) iunit 
  close(iunit)

  c = c_tmp
  Occ = Occ_tmp

  deallocate(Occ_tmp,c_tmp)

end subroutine 
