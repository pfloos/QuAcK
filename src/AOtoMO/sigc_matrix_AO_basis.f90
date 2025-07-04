subroutine sigc_AO_basis(nBas,nOrb,c,U_QP,eqsGW_state,ERI,Sigc)

! Compute Sigma_c matrix in the AO basis

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  double precision,intent(in)   :: U_QP(nOrb*2,nOrb*2)
  double precision,intent(in)   :: eqsGW_state(nOrb*2)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  integer                       :: iorb
  integer                       :: ifreq
  integer                       :: nOrb2
  integer                       :: nfreq=10

  double precision,allocatable  :: R(:,:) 
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  complex*16,allocatable        :: Chi0_hehe(:,:)

! Output variables

  double precision,intent(out)  :: Sigc(nBas,nBas)

  nOrb2 = 2*nOrb
  Sigc(:,:) = 0d0

! Building ERIs in MO basis

  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI,ERI_MO)

! Building Chi0_hehe and Chi0_hhee

  allocate(Chi0_hehe(nOrb*nOrb,nOrb*nOrb))
  do ifreq=1,nfreq
    Chi0_hehe(:,:) = (0d0,0d0)

    

  enddo
 
  

  deallocate(Chi0_hehe)

write(*,'(*(f10.5))') eqsGW_state(:)
write(*,*) ' ' 
do iorb=1,nOrb2
 write(*,'(*(f10.5))') U_QP(iorb,:)
enddo

write(*,*) 'R'
allocate(R(nOrb2,nOrb2))
R(:,:)     = 0d0
do iorb=1,nOrb
 R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
enddo
do iorb=1,nOrb2
  write(*,'(*(f10.5))') R(iorb,:)
enddo
deallocate(R)

  deallocate(ERI_MO)

end subroutine 

! ---

