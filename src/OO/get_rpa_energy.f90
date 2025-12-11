subroutine get_rpa_energy(O,V,N,nS,Hc,c,ERI_AO,ERPA)
      
! Compute RPA energy for given orbitals c

  implicit none
  include 'parameters.h'
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Hc(N,N),c(N,N)
  double precision,intent(in)   :: ERI_AO(N,N,N,N)

! Local variables
  double precision,allocatable  :: rdm1_hf(:,:)
  double precision,allocatable  :: rdm2_hf(:,:,:,:)
  double precision,allocatable  :: PHF(:,:),J(:,:),FHF(:,:),K(:,:),h(:,:),F(:,:)
  double precision,allocatable  :: X(:,:),Y(:,:),XpY(:,:),XmY(:,:),Om(:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: Aph(:,:),Bph(:,:)
  double precision              :: EcRPA
  logical                       :: dRPA = .true.
  integer                       :: isp_W = 1
  logical                       :: TDA_W = .false.


! Output variables

  double precision,intent(out)  :: ERPA

! Allocate stuff
allocate(rdm1_hf(N,N), rdm2_hf(N,N,N,N), PHF(N,N), J(N,N), K(N,N), FHF(N,N), &
         X(nS,nS), Y(nS,nS), XpY(nS,nS), XmY(nS,nS), Om(nS), &
         ERI_MO(N,N,N,N), h(N,N), F(N,N), Aph(nS,nS), Bph(nS,nS))

  ! Compute Fock operator
  PHF(:,:) = 2d0 * matmul(c(:,1:O), transpose(c(:,1:O))) 
  J(:,:) = 0d0
  K(:,:) = 0d0
  call Hartree_matrix_AO_basis(N,PHF,ERI_AO,J)
  call exchange_matrix_AO_basis(N,PHF,ERI_AO,K)
  FHF(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)
  call AOtoMO(N,N,c,FHF,F)
  call AOtoMO(N,N,c,Hc,h)
  call AOtoMO_ERI_RHF(N,N,c,ERI_AO,ERI_MO)

  ! Compute screening
  call OO_phRLR_A(isp_W,dRPA,N,0,O,V,0,nS,1d0,F,ERI_MO,Aph)
  call phRLR_B(isp_W,dRPA,N,0,O,V,0,nS,1d0,ERI_MO,Bph)
  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
  X = transpose(0.5*(XpY + XmY))
  Y = transpose(0.5*(XpY - XmY))

  ! Compute energy
  call RG0W0_rdm2_hf(O,V,N,nS,rdm2_hf)
  call RG0W0_rdm1_hf(O,V,N,nS,rdm1_hf)
  call energy_from_rdm(N,h,ERI_MO,rdm1_hf,rdm2_hf,ERPA,.false.)
  ERPA = ERPA + EcRPA

deallocate(rdm1_hf,rdm2_hf,PHF,J,K,X,Y,XpY,XmY,Om,ERI_MO,h,F,Aph,Bph,FHF)
end subroutine
