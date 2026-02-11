subroutine complex_phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response for complex unrestricted formalism using the c-product

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  complex*16,intent(in)         :: Aph(nSt,nSt)
  complex*16,intent(in)         :: Bph(nSt,nSt)
  
! Local variables

  complex*16,external           :: complex_trace_matrix
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: OmOmminus(:)

! Output variables

  complex*16,intent(out)        :: EcRPA
  complex*16,intent(out)        :: Om(nSt)
  complex*16,intent(out)        :: XpY(nSt,nSt)
  complex*16,intent(out)        :: XmY(nSt,nSt)



! Tamm-Dancoff approximation

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call complex_diagonalize_matrix(nSt,XpY,Om)
    call complex_orthogonalize_matrix(nSt,XpY)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    allocate(RPA_matrix(2*nSt,2*nSt),OmOmminus(2*nSt))
    RPA_matrix(1:nSt,1:nSt) = Aph(:,:)
    RPA_matrix(1:nSt,nSt+1:2*nSt) = Bph(:,:)
    RPA_matrix(nSt+1:2*nSt,1:nSt) = -Bph(:,:)
    RPA_matrix(nSt+1:2*nSt,nSt+1:2*nSt) = -Aph(:,:)
    call complex_diagonalize_matrix_without_sort(2*nSt,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nSt,OmOmminus,RPA_matrix)
    call complex_normalize_RPA(nSt,RPA_matrix)
    Om(:) = OmOmminus(1:nSt)
    if(maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt))) > 1e-12) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt)))
    end if
    if(minval(real(Om(:))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    XpY(:,:) = transpose(RPA_matrix(1:nSt,1:nSt) + RPA_matrix(nSt+1:2*nSt,1:nSt)) 
    XmY(:,:) = transpose(RPA_matrix(1:nSt,1:nSt) - RPA_matrix(nSt+1:2*nSt,1:nSt))
    print *, "Om+"
    call complex_vecout(nSt,OmOmminus)
    print *, "Om-"
    call complex_vecout(nSt,OmOmminus(nSt+1:2*nSt))
    deallocate(RPA_matrix,OmOmminus) 

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - complex_trace_matrix(nSt,Aph))

end subroutine 
