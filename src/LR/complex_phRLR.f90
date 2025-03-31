subroutine phRLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: Aph(nS,nS)
  complex*16,intent(in)         :: Bph(nS,nS)

  ! Local variables

  complex*16                    :: complex_trace_matrix
  complex*16                    :: t1, t2
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: Z(:,:)
  complex*16,allocatable        :: tmp(:,:)
  complex*16,allocatable        :: OmOmminus(:)

! Output variables

  complex*16,intent(out)        :: EcRPA
  complex*16,intent(out)        :: Om(nS)
  complex*16,intent(out)        :: XpY(nS,nS)
  complex*16,intent(out)        :: XmY(nS,nS)



! Tamm-Dancoff approximation

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call complex_diagonalize_matrix(nS,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    allocate(RPA_matrix(2*nS,2*nS),OmOmminus(2*nS))

    call complex_diagonalize_matrix(2*nS,RPA_matrix,Om)
    Om(:) = OmOmminus(1:nS)
    call complex_vecout(nS,Om)
    call complex_vecout(nS,Om(nS+1:2*nS))
    if (maxval(abs(Om(:)+Om(nS+1:2*nS))) > 1e10) &
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
    if(minval(abs(Om(:))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    XpY(:,:) = RPA_matrix(1:nS,1:nS) + RPA_matrix(nS+1:2*nS,nS+1:2*nS) 
    XmY(:,:) = RPA_matrix(1:nS,1:nS) - RPA_matrix(nS+1:2*nS,nS+1:2*nS) 

    deallocate(RPA_matrix)
 
  end if

  ! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - complex_trace_matrix(nS,Aph))

end subroutine 
