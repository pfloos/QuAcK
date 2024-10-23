subroutine phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response for unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: Aph(nSt,nSt)
  double precision,intent(in)   :: Bph(nSt,nSt)
  
! Local variables

  double precision,external     :: trace_matrix
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

  allocate(ApB(nSt,nSt),AmB(nSt,nSt),AmBSq(nSt,nSt),AmBIv(nSt,nSt),Z(nSt,nSt))

! Tamm-Dancoff approximation

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nSt,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    ApB(:,:) = Aph(:,:) + Bph(:,:)
    AmB(:,:) = Aph(:,:) - Bph(:,:)

  ! Diagonalize linear response matrix

    call diagonalize_matrix(nSt,AmB,Om)

    if(minval(Om) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')

    call ADAt(nSt,AmB,1d0*sqrt(Om),AmBSq)
    call ADAt(nSt,AmB,1d0/sqrt(Om),AmBIv)
 
    Z = matmul(AmBSq,matmul(ApB,AmBSq))
 
    call diagonalize_matrix(nSt,Z,Om)

    if(minval(Om) < 0d0) & 
      call print_warning('You may have instabilities in linear response: negative excitations!!')
    
    Om = sqrt(Om)
 
    XpY = matmul(transpose(Z),AmBSq)
    call DA(nSt,1d0/sqrt(Om),XpY)
 
    XmY = matmul(transpose(Z),AmBIv)
    call DA(nSt,1d0*sqrt(Om),XmY)

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

end subroutine 
