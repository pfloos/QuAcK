subroutine CVS_phULR(TDA,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response for complex unrestricted formalism using the c-product

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
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: OmOmminus(:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Om(nSt)
  double precision,intent(out)  :: XpY(nSt,nSt)
  double precision,intent(out)  :: XmY(nSt,nSt)

! Tamm-Dancoff approximation

  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nSt,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

  else

    allocate(RPA_matrix(2*nSt,2*nSt),OmOmminus(2*nSt))
    
    RPA_matrix(1:nSt,1:nSt)             = cmplx( Aph(:,:),0d0,kind=8)
    RPA_matrix(1:nSt,nSt+1:2*nSt)       = cmplx( Bph(:,:),0d0,kind=8)
    RPA_matrix(nSt+1:2*nSt,1:nSt)       = cmplx(-Bph(:,:),0d0,kind=8)
    RPA_matrix(nSt+1:2*nSt,nSt+1:2*nSt) = cmplx(-Aph(:,:),0d0,kind=8)
    
    call complex_diagonalize_matrix_without_sort(2*nSt,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nSt,OmOmminus,RPA_matrix)
    call complex_normalize_RPA(nSt,RPA_matrix)
    call rotate_vectors_to_real_axis(2*nSt,nSt,RPA_matrix(1:2*nSt,1:nSt))

    Om(:) = real(OmOmminus(1:nSt))
    if(maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt))) > 1e-8) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt)))
    end if
    if(minval(real(OmOmminus(1:nSt))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    if(maxval(abs(aimag(OmOmminus(:))))>1d-8)&
      call print_warning('You may have instabilities in linear response: complex excitation eigenvalue !')
    if(maxval(aimag(RPA_matrix(1:2*nSt,1:nSt)))>1d-8) then
      call print_warning('You may have instabilities in linear response: complex transition vectors !')
      print *, "Max imag value X+Y:",maxval(aimag(transpose(RPA_matrix(1:nSt,1:nSt) + RPA_matrix(nSt+1:2*nSt,1:nSt))))
      print *, "Max imag value X-Y:",maxval(aimag(transpose(RPA_matrix(1:nSt,1:nSt) - RPA_matrix(nSt+1:2*nSt,1:nSt))))
    end if
    XpY(:,:) = transpose(real(RPA_matrix(1:nSt,1:nSt) + RPA_matrix(nSt+1:2*nSt,1:nSt))) 
    XmY(:,:) = transpose(real(RPA_matrix(1:nSt,1:nSt) - RPA_matrix(nSt+1:2*nSt,1:nSt))) 
    deallocate(RPA_matrix,OmOmminus) 

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

end subroutine 
