subroutine unrestricted_linear_response(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda, & 
                                        e,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,Omega,XpY,XmY)

! Compute linear response for unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  double precision,intent(in)   :: eta
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
  double precision,allocatable  :: A(:,:)
  double precision,allocatable  :: B(:,:)
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Omega(nSt)
  double precision,intent(out)  :: XpY(nSt,nSt)
  double precision,intent(out)  :: XmY(nSt,nSt)

! Memory allocation

  allocate(A(nSt,nSt),B(nSt,nSt),ApB(nSt,nSt),AmB(nSt,nSt),AmBSq(nSt,nSt),AmBIv(nSt,nSt),Z(nSt,nSt))

! Build A and B matrices 

  call unrestricted_linear_response_A_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda,e, & 
                                             ERI_aaaa,ERI_aabb,ERI_bbbb,A)

  if(BSE) & 
    call unrestricted_Bethe_Salpeter_A_matrix(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda,e, & 
                                              ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,A)

! Tamm-Dancoff approximation

  if(TDA) then

    B(:,:)   = 0d0
    XpY(:,:) = A(:,:)
    call diagonalize_matrix(nSt,XpY,Omega)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    call unrestricted_linear_response_B_matrix(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,lambda, & 
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,B)

    if(BSE) &
      call unrestricted_Bethe_Salpeter_B_matrix(ispin,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,nS_sc,lambda, & 
                                                ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,B)

  ! Build A + B and A - B matrices 

    ApB = A + B
    AmB = A - B

  ! Diagonalize linear response matrix

   call diagonalize_matrix(nSt,AmB,Omega)

    if(minval(Omega) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')

  ! do ia=1,nSt
  !   if(Omega(ia) < 0d0) Omega(ia) = 0d0
  ! end do

    call ADAt(nSt,AmB,1d0*sqrt(Omega),AmBSq)
    call ADAt(nSt,AmB,1d0/sqrt(Omega),AmBIv)
 
    Z = matmul(AmBSq,matmul(ApB,AmBSq))
 
    call diagonalize_matrix(nSt,Z,Omega)

    if(minval(Omega) < 0d0) & 
      call print_warning('You may have instabilities in linear response: negative excitations!!')
    
  ! do ia=1,nSt
  !   if(Omega(ia) < 0d0) Omega(ia) = 0d0
  ! end do

    Omega = sqrt(Omega)
 
    XpY = matmul(transpose(Z),AmBSq)
    call DA(nSt,1d0/sqrt(Omega),XpY)
 
    XmY = matmul(transpose(Z),AmBIv)
    call DA(nSt,1d0*sqrt(Omega),XmY)

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Omega) - trace_matrix(nSt,A))

end subroutine unrestricted_linear_response
