subroutine RGWC(dotest,nBas,eGW,Z)

! Perform GW+C calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: nBas

  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Z(nBas)

! Local variables

  integer                       :: p

  double precision,allocatable  :: eGWC(:)
  double precision,allocatable  :: ZC(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted GW+C Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Memory allocation

  allocate(eGWC(nBas),ZC(nBas))

! GW+C weights

  ZC(:) = exp(1d0 - 1d0/Z(:))

! GW+C quasiparticles

  eGWC(:) = eGW(:)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)' GW+C calculation '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_GW (eV)','|','e_GW+C (eV)','|','Z','|','ZC','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eGW(p)*HaToeV,'|',eGWC(p)*HaToeV,'|',Z(p),'|',ZC(p),'|'
  end do

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Testing zone

! if(dotest) then

!   call dump_test_value('R','G0W0 correlation energy',EcRPA)
!   call dump_test_value('R','G0W0 HOMO energy',eGW(nO))
!   call dump_test_value('R','G0W0 LUMO energy',eGW(nO+1))

! end if

end subroutine 
