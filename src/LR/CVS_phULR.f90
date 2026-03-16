subroutine CVS_phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

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
  double precision              :: eF
  integer                       :: i

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Om(nSt)
  double precision,intent(out)  :: XpY(nSt,nSt)
  double precision,intent(out)  :: XmY(nSt,nSt)

! Memory allocation

  allocate(ApB(nSt,nSt),AmB(nSt,nSt))

! Tamm-Dancoff approximation

!  eF       = 20d0
!  do i=1,nSt
!    Aph(i,i) = Aph(i,i) + 2*eF
!  enddo

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nSt,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    ApB(:,:) = Aph(:,:) + Bph(:,:) 
    AmB(:,:) = Aph(:,:) - Bph(:,:)


  ! Diagonalize linear response matrix
   
!    ApB = matmul(ApB,AmB)
!    ApB = ApB - 4*eF*Aph
!    call diagonalize_general_matrix(nSt,ApB,Om,AmB)
!    call vecout(nSt,sqrt(Om + 4*eF**2))
    ApB = matmul(ApB,AmB)
    call diagonalize_general_matrix(nSt,ApB,Om,XmY)
    Om = sqrt(Om)
    XpY = matmul(AmB,XmY)
    XpY = transpose(XpY)
    call DA(nSt,1d0/Om,XpY)
    XpY = transpose(XpY)
    
    ApB = matmul(transpose(XmY),XpY)
    do i=1,nSt
      if(ApB(i,i)<0d0) then
        Om(i)    = - Om(i)
        XmY(:,i) = - XmY(:,i)
      end if
    enddo
    call sort_eigenvalues_vec_vec(nSt,Om,XmY,XpY)
    ApB = matmul(transpose(XmY),XpY)
    call orthogonalize_matrix(1,nSt,ApB,AmB)
    XmY = matmul(XmY,AmB)
    XpY = matmul(XpY,AmB)
    XmY = transpose(XmY) 
    XpY = transpose(XpY) 
    
  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

end subroutine 
