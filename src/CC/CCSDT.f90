subroutine CCSDT(nC,nO,nV,nR,eO,eV,OOVV,VVVO,VOOO,t1,t2,EcCCT)

! Compute the (T) correction of the CCSD(T) energy

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR

  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: eV(nV)

  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: VVVO(nV,nV,nV,nO)
  double precision,intent(in)   :: VOOO(nV,nO,nO,nO)

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

! Local variables

  double precision,allocatable  :: delta_OOOVVV(:,:,:,:,:,:)
  double precision,allocatable  :: ub(:,:,:,:,:,:)
  double precision,allocatable  :: ubb(:,:,:,:,:,:)

! Output variables

  double precision,intent(out)  :: EcCCT

! Memory allocation

  allocate(delta_OOOVVV(nO,nO,nO,nV,nV,nV),ub(nO,nO,nO,nV,nV,nV),ubb(nO,nO,nO,nV,nV,nV))

! Form CCSD(T) quantities

  call form_delta_OOOVVV(nC,nO,nV,nR,eO,eV,delta_OOOVVV) 

  call form_ub(nC,nO,nV,nR,OOVV,t1,ub)

  call form_ubb(nC,nO,nV,nR,VVVO,VOOO,t2,ubb)

  call form_T(nC,nO,nV,nR,delta_OOOVVV,ub,ubb,EcCCT)

end subroutine 
