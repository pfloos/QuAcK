subroutine complex_phRLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response

  implicit none
  include 'parameters.h'

  ! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nS
  complex*16,intent(inout)      :: Aph(nS,nS)
  complex*16,intent(in)         :: Bph(nS,nS)

  ! Local variables

  complex*16,external           :: complex_trace_matrix
  complex*16                    :: t1, t2
  complex*16,allocatable        :: RPA_matrix(:,:)
  complex*16,allocatable        :: OmOmminus(:)

  double precision,allocatable  :: RPA_matrix_real(:,:)
  double precision,allocatable  :: U(:,:)
  double precision,allocatable  :: Vt(:,:)
  double precision,allocatable  :: D(:)
  double precision              :: s
  logical                       :: reg = .true.

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
    XpY = transpose(XpY)
    XmY(:,:) = XpY(:,:)

  else

    allocate(RPA_matrix(2*nS,2*nS),OmOmminus(2*nS))
    RPA_matrix(1:nS,1:nS) = Aph(:,:)
    RPA_matrix(1:nS,nS+1:2*nS) = Bph(:,:)
    RPA_matrix(nS+1:2*nS,1:nS) = -Bph(:,:)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)
    call complex_matout(2*nS,2*nS,RPA_matrix)

    if(reg) then
      print *, "SVD Regularisation of RPA problem"
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Singularvalue decomp  !
      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! SVD of RPA matrix
      allocate(RPA_matrix_real(2*nS,2*nS),U(2*nS,2*nS),Vt(2*nS,2*nS),D(2*nS))
      s = 10d0
      RPA_matrix_real = real(RPA_matrix)
      call svd(2*nS,RPA_matrix_real,U,D,Vt)
      print*, "Singular values of RPA matrix"
      call vecout(2*nS,D)
      ! Approx rpa matrix
      D(:) = D(:)*(1-exp(-D(:)**2*s))
      call AD(2*nS,U,D)
      RPA_matrix_real = matmul(U,Vt)
      RPA_matrix = cmplx(RPA_matrix_real,0d0,kind=8)
      deallocate(RPA_matrix_real,U,Vt,D)
      !! SVD of TDA problem
      !allocate(RPA_matrix_real(nS,nS),U(nS,nS),Vt(nS,nS),D(nS))
      allocate(RPA_matrix_real(2*nS,2*nS),U(2*nS,2*nS),Vt(2*nS,2*nS),D(2*nS))
      RPA_matrix_real = 0d0
      RPA_matrix_real(1:nS,1:nS) = real(Aph)
      RPA_matrix_real(nS+1:2*nS,nS+1:2*nS) = -real(Aph)
      call svd(2*nS,RPA_matrix_real,U,D,Vt)
      print*, "Singular values of Aph matrix"
      call vecout(2*nS,D)
      ! Approx rpa matrix
      D(:) = D(:)*(1-exp(-D(:)**2*s))
      call AD(2*nS,U,D)
      RPA_matrix_real = matmul(U,Vt)
      Aph             = cmplx(RPA_matrix_real,0d0,kind=8)
      deallocate(RPA_matrix_real,U,Vt,D)
   end if 


    !!!!!!!!!!!!!!!!!!!!!!!!!!

    call complex_diagonalize_matrix_without_sort(2*nS,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nS,OmOmminus,RPA_matrix)
    call complex_normalize_RPA(nS,RPA_matrix)
    call complex_vecout(2*nS,OmOmminus)
    Om(:) = OmOmminus(1:nS)
    if(maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS))) > 1e-12) then
      call print_warning('We dont find a Om and -Om structure as solution of the RPA. There might be a problem somewhere.')
      write(*,*) "Maximal difference :", maxval(abs(OmOmminus(1:nS)+OmOmminus(nS+1:2*nS)))
      call complex_vecout(nS,OmOmminus(1:nS)+OmOmminus(nS+1:2*nS))
    end if
    if(minval(real(Om(:))) < 0d0) &
      call print_warning('You may have instabilities in linear response: A-B is not positive definite!!')
    XpY(:,:) = transpose(RPA_matrix(1:nS,1:nS) + RPA_matrix(nS+1:2*nS,1:nS)) 
    XmY(:,:) = transpose(RPA_matrix(1:nS,1:nS) - RPA_matrix(nS+1:2*nS,1:nS))
    deallocate(RPA_matrix,OmOmminus) 
  end if

  ! Compute the RPA correlation energy

  EcRPA = (0.5d0,0.d0)*(sum(Om) - complex_trace_matrix(nS,Aph))

end subroutine 
