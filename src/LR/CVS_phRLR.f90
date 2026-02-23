subroutine CVS_phRLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response for complex unrestricted formalism using the c-product

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nS
  double precision,intent(in)   :: Aph(nS,nS)
  double precision,intent(in)   :: Bph(nS,nS)
  
! Local variables

  double precision,external     :: trace_matrix
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: OmOmminus(:)

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

    allocate(RPA_matrix(2*nS,2*nS),OmOmminus(2*nS))
    RPA_matrix(1:nS,1:nS) = cmplx(Aph(:,:),0d0,kind=8)
    RPA_matrix(1:nS,nS+1:2*nS) = cmplx(Bph(:,:),0d0,kind=8)
    RPA_matrix(nS+1:2*nS,1:nS) = cmplx(-Bph(:,:),0d0,kind=8)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = cmplx(-Aph(:,:),0d0,kind=8)
    call complex_diagonalize_matrix_without_sort(2*nS,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nS,OmOmminus,RPA_matrix)
    call complex_normalize_RPA(nS,RPA_matrix)

    Om(:) = real(OmOmminus(1:nS))
    if(maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS))) > 1e-12) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS)))
    end if
    if(minval(real(OmOmminus(:))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    if(maxval(abs(aimag(OmOmminus(:))))>1d-8)&
      call print_warning('You may have instabilities in linear response: complex excitation eigenvalue !')
    XpY(:,:) = transpose(real(RPA_matrix(1:nS,1:nS) + RPA_matrix(nS+1:2*nS,1:nS))) 
    XmY(:,:) = transpose(real(RPA_matrix(1:nS,1:nS) - RPA_matrix(nS+1:2*nS,1:nS)))
    
    deallocate(RPA_matrix,OmOmminus) 

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nS,Aph))

end subroutine 
