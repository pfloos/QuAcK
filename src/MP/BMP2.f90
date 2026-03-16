subroutine BMP2(nBas,nOrb,cHFB,Hc,S,ERI,chem_pot,sigma,U_QP,ERHFB,EcMP2)

! Perform second-order Bogoliubov Moller-Plesset calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  double precision,intent(in)   :: sigma,chem_pot
  double precision,intent(in)   :: ERHFB
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: S(nBas,nBas) 
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: U_QP(2*nOrb,2*nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: iorb,jorb,korb,lorb
  integer                       :: nOrb2,nOrb3,nOrb4

  double precision              :: trace_1rdm
  double precision              :: EcdMP2,EcxMP2

  double precision,allocatable  :: eQP_sw(:)
  double precision,allocatable  :: F(:,:),J(:,:),K(:,:),Delta(:,:)
  double precision,allocatable  :: P(:,:),Panom(:,:),R(:,:)
  double precision,allocatable  :: U_QP_sw(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: ERI_MO_sw(:,:,:,:)
  double precision,allocatable  :: Ua(:,:),Va(:,:),U(:,:),V(:,:)

! Output variables

  double precision,intent(out)  :: EcMP2

! Hello world

  write(*,*)
  write(*,*)'******************************'
  write(*,*)'* Bogoliubov MP2 Calculation *'
  write(*,*)'******************************'
  write(*,*)

  nOrb2=2*nOrb
  nOrb3=nOrb2+nOrb
  nOrb4=2*nOrb2
  allocate(U_QP_sw(nOrb4,nOrb4))
  allocate(eQP_sw(nOrb4))
  allocate(ERI_MO_sw(nOrb2,nOrb2,nOrb2,nOrb2))

! Build P and Panom -> U_QP_sw, eQP_sw
  ! Allocate arrays
  allocate(P(nBas,nBas))
  allocate(Panom(nBas,nBas))
  allocate(R(nOrb2,nOrb2))
  allocate(F(nBas,nBas))
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Delta(nBas,nBas))
  trace_1rdm=0d0
  R=0d0
  do iorb=1,nOrb
   R(:,:) = R(:,:) + matmul(U_QP(:,iorb:iorb),transpose(U_QP(:,iorb:iorb))) 
  enddo
  do iorb=1,nOrb
   trace_1rdm = trace_1rdm + R(iorb,iorb) 
  enddo
  P=2d0*matmul(matmul(cHFB,R(1:nOrb,1:nOrb)),transpose(cHFB))
  Panom=matmul(matmul(cHFB,R(1:nOrb,nOrb+1:nOrb2)),transpose(cHFB))
  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)
  call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) - chem_pot*S(:,:)
  U_QP_sw=0d0
  ! 1 row-block
  U_QP_sw(1:nOrb      ,1:nOrb      )              = matmul(transpose(cHFB),matmul(F,cHFB))
  U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4) = matmul(transpose(cHFB),matmul(Delta,cHFB))
  ! 2 row-block
  U_QP_sw(nOrb+1:nOrb2,nOrb+1:nOrb2)        =  U_QP_sw(1:nOrb,1:nOrb)
  U_QP_sw(nOrb+1:nOrb2,nOrb2+1:nOrb3) = -U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4)
  ! 3 row-block
  U_QP_sw(nOrb2+1:nOrb3,nOrb2+1:nOrb3)   = -U_QP_sw(1:nOrb,1:nOrb)
  U_QP_sw(nOrb2+1:nOrb3,nOrb+1:nOrb2)          =  U_QP_sw(nOrb+1:nOrb2,nOrb2+1:nOrb3)
  ! 4 row-block 
  U_QP_sw(nOrb3+1:nOrb4,nOrb3+1:nOrb4) = -U_QP_sw(1:nOrb,1:nOrb)
  U_QP_sw(nOrb3+1:nOrb4,1:nOrb                   ) =  U_QP_sw(1:nOrb      ,nOrb3+1:nOrb4)
  call diagonalize_matrix(nOrb4,U_QP_sw,eQP_sw)
  !  ! Test R is in NO and spin-with
  !  deallocate(R)
  !  allocate(R(nOrb4,nOrb4))
  !  trace_1rdm=0d0
  !  R=0d0
  !  do iorb=1,nOrb2
  !   R(:,:) = R(:,:) + matmul(U_QP_sw(:,iorb:iorb),transpose(U_QP_sw(:,iorb:iorb))) 
  !  enddo
  !  do iorb=1,nOrb2
  !   trace_1rdm = trace_1rdm + R(iorb,iorb) 
  !  enddo
  !  write(*,'(a,f10.5)') ' Trace 1D ',trace_1rdm
  !  write(*,*) ' spin-with eQP'
  !  do iorb=1,nOrb4
  !   write(*,'(*(f10.5))') eQP_sw(iorb)
  !  enddo
  !  write(*,*) ' R in spin-with and NO basis'
  !  do iorb=1,nOrb4
  !   write(*,'(*(f10.5))') R(iorb,:)
  !  enddo
  deallocate(P,Panom,R)
  deallocate(F,J,K,Delta)

! Allocate NO spin-with ERIs
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,cHFB,ERI,ERI_MO)
  ERI_MO_sw=0d0
  ERI_MO_sw(      1:nOrb,      1:nOrb,      1:nOrb,      1:nOrb)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! a a a a 
  ERI_MO_sw(nOrb+1:nOrb2,nOrb+1:nOrb2,nOrb+1:nOrb2,nOrb+1:nOrb2)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! b b b b 
  ERI_MO_sw(      1:nOrb,nOrb+1:nOrb2,      1:nOrb,nOrb+1:nOrb2)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! a b a b 
  ERI_MO_sw(nOrb+1:nOrb2,      1:nOrb,nOrb+1:nOrb2,      1:nOrb)=ERI_MO(1:nOrb,1:nOrb,1:nOrb,1:nOrb) ! b a b a 
  deallocate(ERI_MO)

  deallocate(ERI_MO_sw)
  deallocate(U_QP_sw)
  deallocate(eQP_sw)

end subroutine 
