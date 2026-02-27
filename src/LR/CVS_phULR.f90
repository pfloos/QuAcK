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
  double precision,allocatable  :: RPA_matrix(:,:)
  complex*16,allocatable        :: complex_RPA_matrix(:,:)
  double precision,allocatable  :: vectors(:,:)
  double precision,allocatable  :: OmOmminus(:)
  integer                       :: i

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

    allocate(RPA_matrix(2*nSt,2*nSt),vectors(2*nSt,2*nSt),OmOmminus(2*nSt))
    
    RPA_matrix(1:nSt,1:nSt)             =  Aph(:,:)
    RPA_matrix(1:nSt,nSt+1:2*nSt)       =  Bph(:,:)
    RPA_matrix(nSt+1:2*nSt,1:nSt)       = -Bph(:,:)
    RPA_matrix(nSt+1:2*nSt,nSt+1:2*nSt) = -Aph(:,:)
    
    call diagonalize_general_matrix(2*nSt,RPA_matrix,OmOmminus,vectors)
    
    RPA_matrix(:,:) = vectors(:,:)
    
    deallocate(vectors)
    
    call sort_eigenvalues_RPA(2*nSt,OmOmminus,RPA_matrix)
    
    allocate(complex_RPA_matrix(2*nSt,2*nSt))
    
    complex_RPA_matrix = cmplx(0.0d0, 0.0d0, kind=8)
    complex_RPA_matrix = cmplx(RPA_matrix(:,:),0d0,kind=8) 
    
    deallocate(RPA_matrix)
    
    call complex_normalize_RPA(nSt,complex_RPA_matrix)
    
    if(maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt))) > 1e-8) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nSt)+OmOmminus(nSt+1:2*nSt)))
    end if

    Om(:) = OmOmminus(1:nSt)
    
    deallocate(OmOmminus)
    
    ! Set transition vectors of zero modes to 0
    do i=1,nSt
      if(abs(Om(i))<1d-8) then
         complex_RPA_matrix(:,i) = (0d0,0d0)
         complex_RPA_matrix(:,i+nSt) = (0d0,0d0)
      end if
    end do


    if(maxval(aimag(complex_RPA_matrix(1:2*nSt,1:nSt)))>1d-8) then
      call print_warning('You may have instabilities in linear response: complex transition vectors !')
      print *, "Max imag value X+Y:",maxval(aimag(transpose(complex_RPA_matrix(1:nSt,1:nSt) + complex_RPA_matrix(nSt+1:2*nSt,1:nSt)))),&
               "at"                 ,maxloc(aimag(transpose(complex_RPA_matrix(1:nSt,1:nSt) + complex_RPA_matrix(nSt+1:2*nSt,1:nSt))))
 
      print *, "Max imag value X-Y:",maxval(aimag(transpose(complex_RPA_matrix(1:nSt,1:nSt) - complex_RPA_matrix(nSt+1:2*nSt,1:nSt)))),&
               "at"                 ,maxloc(aimag(transpose(complex_RPA_matrix(1:nSt,1:nSt) - complex_RPA_matrix(nSt+1:2*nSt,1:nSt))))
    end if

    XpY(:,:) = transpose(real(complex_RPA_matrix(1:nSt,1:nSt) + complex_RPA_matrix(nSt+1:2*nSt,1:nSt))) 
    XmY(:,:) = transpose(real(complex_RPA_matrix(1:nSt,1:nSt) - complex_RPA_matrix(nSt+1:2*nSt,1:nSt))) 
    
    deallocate(complex_RPA_matrix) 
  
  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

end subroutine 
