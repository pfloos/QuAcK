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
  complex*16,allocatable        :: complex_RPA_matrix(:,:)
  double precision,allocatable  :: RPA_matrix(:,:)
  double precision,allocatable  :: vectors(:,:)
  double precision,allocatable  :: OmOmminus(:)
  integer                       :: i

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

    allocate(RPA_matrix(2*nS,2*nS),vectors(2*nS,2*nS),OmOmminus(2*nS))
    
    RPA_matrix(1:nS,1:nS)             =  Aph(:,:)
    RPA_matrix(1:nS,nS+1:2*nS)       =  Bph(:,:)
    RPA_matrix(nS+1:2*nS,1:nS)       = -Bph(:,:)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)
    
    call diagonalize_general_matrix(2*nS,RPA_matrix,OmOmminus,vectors)
    
    RPA_matrix(:,:) = vectors(:,:)
    
    deallocate(vectors)
    
    call sort_eigenvalues_RPA(2*nS,OmOmminus,RPA_matrix)
    
    allocate(complex_RPA_matrix(2*nS,2*nS))
    
    complex_RPA_matrix = cmplx(0.0d0, 0.0d0, kind=8)
    complex_RPA_matrix = cmplx(RPA_matrix(:,:),0d0,kind=8) 
    
    deallocate(RPA_matrix)
    
    call complex_normalize_RPA(nS,complex_RPA_matrix)
    
    if(maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS))) > 1e-8) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS)))
    end if

    Om(:) = OmOmminus(1:nS)
    
    deallocate(OmOmminus)
    
    ! Set transition vectors of zero modes to 0
    do i=1,nS
      if(abs(Om(i))<1d-8) then
         complex_RPA_matrix(:,i) = (0d0,0d0)
         complex_RPA_matrix(:,i+nS) = (0d0,0d0)
      end if
    end do


    if(maxval(aimag(complex_RPA_matrix(1:2*nS,1:nS)))>1d-8) then
      call print_warning('You may have instabilities in linear response: complex transition vectors !')
      print *, "Max imag value X+Y:",maxval(aimag(transpose(complex_RPA_matrix(1:nS,1:nS) + complex_RPA_matrix(nS+1:2*nS,1:nS)))),&
               "at"                 ,maxloc(aimag(transpose(complex_RPA_matrix(1:nS,1:nS) + complex_RPA_matrix(nS+1:2*nS,1:nS))))
 
      print *, "Max imag value X-Y:",maxval(aimag(transpose(complex_RPA_matrix(1:nS,1:nS) - complex_RPA_matrix(nS+1:2*nS,1:nS)))),&
               "at"                 ,maxloc(aimag(transpose(complex_RPA_matrix(1:nS,1:nS) - complex_RPA_matrix(nS+1:2*nS,1:nS))))
    end if

    XpY(:,:) = transpose(real(complex_RPA_matrix(1:nS,1:nS) + complex_RPA_matrix(nS+1:2*nS,1:nS))) 
    XmY(:,:) = transpose(real(complex_RPA_matrix(1:nS,1:nS) - complex_RPA_matrix(nS+1:2*nS,1:nS))) 
    
    deallocate(complex_RPA_matrix) 

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*(sum(Om) - trace_matrix(nS,Aph))

end subroutine 
