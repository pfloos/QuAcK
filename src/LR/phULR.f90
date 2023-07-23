subroutine phULR(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,e, & 
                 ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,Om,XpY,XmY)

! Compute linear response for unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: nS_sc
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

  double precision,intent(in)   :: OmRPA(nS_sc)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS_sc,nspin)
  
! Local variables

  double precision,external     :: trace_matrix
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Om(nSt)
  double precision,intent(out)  :: XpY(nSt,nSt)
  double precision,intent(out)  :: XmY(nSt,nSt)

! Memory allocation

  allocate(Aph(nSt,nSt),Bph(nSt,nSt),ApB(nSt,nSt),AmB(nSt,nSt),AmBSq(nSt,nSt),AmBIv(nSt,nSt),Z(nSt,nSt))

! Build A and B matrices 

  call phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda,e,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)

  if(BSE) & 
    call unrestricted_Bethe_Salpeter_A_matrix(ispin,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,e, & 
                                              ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,Aph)

! Tamm-Dancoff approximation

  if(TDA) then

    Bph(:,:)   = 0d0
    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nSt,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    call phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

    if(BSE) &
      call unrestricted_Bethe_Salpeter_B_matrix(ispin,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda, & 
                                                ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,Bph)

  ! Build A + B and A - B matrices 

    ApB(:,:) = Aph(:,:) + Bph(:,:)
    AmB(:,:) = Aph(:,:) - Bph(:,:)

  ! Diagonalize linear response matrix

   call diagonalize_matrix(nSt,AmB,Om)

    if(minval(Om) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')

  ! do ia=1,nSt
  !   if(Om(ia) < 0d0) Om(ia) = 0d0
  ! end do

    call ADAt(nSt,AmB,1d0*sqrt(Om),AmBSq)
    call ADAt(nSt,AmB,1d0/sqrt(Om),AmBIv)
 
    Z = matmul(AmBSq,matmul(ApB,AmBSq))
 
    call diagonalize_matrix(nSt,Z,Om)

    if(minval(Om) < 0d0) & 
      call print_warning('You may have instabilities in linear response: negative excitations!!')
    
  ! do ia=1,nSt
  !   if(Om(ia) < 0d0) Om(ia) = 0d0
  ! end do

    Om = sqrt(Om)
 
    XpY = matmul(transpose(Z),AmBSq)
    call DA(nSt,1d0/sqrt(Om),XpY)
 
    XmY = matmul(transpose(Z),AmBIv)
    call DA(nSt,1d0*sqrt(Om),XmY)

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

end subroutine 
