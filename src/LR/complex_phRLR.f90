subroutine complex_phRLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: Aph(nS,nS)
  complex*16,intent(in)         :: Bph(nS,nS)

  ! Local variables

  complex*16,external           :: complex_trace_matrix
  complex*16                    :: t1, t2
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: Z(:,:)
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
    call complex_orthogonalize_matrix(nS,XpY)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    allocate(RPA_matrix(2*nS,2*nS),OmOmminus(2*nS))
    RPA_matrix(1:nS,1:nS) = Aph(:,:)
    RPA_matrix(1:nS,nS+1:2*nS) = Bph(:,:)
    RPA_matrix(nS+1:2*nS,1:nS) = -Bph(:,:)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)
    call complex_diagonalize_matrix_without_sort(2*nS,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nS,OmOmminus,RPA_matrix)
    call complex_normalize_RPA(nS,RPA_matrix)
    Om(:) = OmOmminus(1:nS)
    if(maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS))) > 1e-12) &
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
    if(minval(abs(Om(:))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    XpY(:,:) = transpose(RPA_matrix(1:nS,1:nS) + RPA_matrix(nS+1:2*nS,1:nS)) 
    XmY(:,:) = transpose(RPA_matrix(1:nS,1:nS) - RPA_matrix(nS+1:2*nS,1:nS))
    deallocate(RPA_matrix,OmOmminus) 
  end if

  ! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - complex_trace_matrix(nS,Aph))

end subroutine 
