subroutine phLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Aph(nS,nS)
  double precision,intent(in)   :: Bph(nS,nS)

  ! Local variables

  double precision              :: trace_matrix
  double precision              :: t1, t2
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: AmBSq(:,:)
  double precision,allocatable  :: AmBIv(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: tmp(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Om(nS)
  double precision,intent(out)  :: XpY(nS,nS)
  double precision,intent(out)  :: XmY(nS,nS)



! Tamm-Dancoff approximation

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nS,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    allocate(ApB(nS,nS),AmB(nS,nS),AmBSq(nS,nS),AmBIv(nS,nS),Z(nS,nS),tmp(nS,nS))

    ApB(:,:) = Aph(:,:) + Bph(:,:)
    AmB(:,:) = Aph(:,:) - Bph(:,:)

  ! Diagonalize linear response matrix

    call diagonalize_matrix(nS,AmB,Om)

    if(minval(Om) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')

    call ADAt(nS,AmB,1d0*dsqrt(Om),AmBSq)
    call ADAt(nS,AmB,1d0/dsqrt(Om),AmBIv)

    call dgemm('N','N',nS,nS,nS,1d0,ApB,size(ApB,1),AmBSq,size(AmBSq,1),0d0,tmp,size(tmp,1))
    call dgemm('N','N',nS,nS,nS,1d0,AmBSq,size(AmBSq,1),tmp,size(tmp,1),0d0,Z,size(Z,1))

!   Z = matmul(AmBSq,matmul(ApB,AmBSq))

    call diagonalize_matrix(nS,Z,Om)

    if(minval(Om) < 0d0) & 
      call print_warning('You may have instabilities in linear response: negative excitations!!')
 
    Om = sqrt(Om)

    call dgemm('T','N',nS,nS,nS,1d0,Z,size(Z,1),AmBSq,size(AmBSq,1),0d0,XpY,size(XpY,1))
    call DA(nS,1d0/dsqrt(Om),XpY)

    call dgemm('T','N',nS,nS,nS,1d0,Z,size(Z,1),AmBIv,size(AmBIv,1),0d0,XmY,size(XmY,1))
    call DA(nS,1d0*dsqrt(Om),XmY)

!   XpY = matmul(transpose(Z),AmBSq)
!   call DA(nS,1d0/sqrt(Om),XpY)

!   XmY = matmul(transpose(Z),AmBIv)
!   call DA(nS,1d0*sqrt(Om),XmY)

    deallocate(ApB,AmB,AmBSq,AmBIv,Z,tmp)
 
  end if

  ! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nS,Aph))

end subroutine 
