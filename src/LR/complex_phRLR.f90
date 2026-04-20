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
  double precision,allocatable  :: TDA_matrix_real(:,:)
  complex*16,allocatable        :: TDA_matrix(:,:),complex_U(:,:),complex_Vt(:,:),X_inv(:,:)
  double precision,allocatable  :: U(:,:)
  double precision,allocatable  :: Vt(:,:)
  double precision,allocatable  :: D(:)
  integer                       :: i 
  double precision              :: thresh_reg = 1.0d0

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
    
    RPA_matrix(1:nS,1:nS)           = Aph(:,:)
    RPA_matrix(1:nS,nS+1:2*nS)      = Bph(:,:)
    RPA_matrix(nS+1:2*nS,1:nS)      = -Bph(:,:)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)

    allocate(RPA_matrix_real(2*nS,2*nS),U(2*nS,2*nS),Vt(2*nS,2*nS),D(2*nS),&
             TDA_matrix_real(2*nS,2*nS),TDA_matrix(2*nS,2*nS))
    allocate(complex_U(2*nS,2*nS),complex_Vt(2*nS,2*nS))
    
    call complex_diagonalize_matrix_without_sort(2*nS,RPA_matrix,OmOmminus)
    call complex_sort_eigenvalues_RPA(2*nS,OmOmminus,RPA_matrix)
    do i=1,nS
      D(i)    =  1d0
      D(i+nS) = -1d0
    enddo
    TDA_matrix = RPA_matrix
    call complex_DA(2*nS,cmplx(D,kind=8),RPA_matrix)
    RPA_matrix = matmul(conjg(transpose(TDA_matrix)),RPA_matrix)
    print  *,"Metrik"
    call complex_matout(2*nS,2*nS,RPA_matrix)
    call complex_svd(2*nS,RPA_matrix,complex_U,D,complex_Vt)
    print*, "SVD of Metrik"
    call vecout(2*nS,D)
    
    ! RPA in other basis
    RPA_matrix(1:nS,1:nS)           = Aph(:,:)
    RPA_matrix(1:nS,nS+1:2*nS)      = Bph(:,:)
    RPA_matrix(nS+1:2*nS,1:nS)      = -Bph(:,:)
    RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)
    print  *,"RPA matrix"
    call complex_matout(2*nS,2*nS,RPA_matrix)
    
    RPA_matrix = matmul(matmul(transpose(conjg(complex_U)),RPA_matrix)&
                                   ,transpose(conjg(complex_Vt)))
    if(abs(D(2*nS))<thresh_reg) then
      RPA_matrix(:,2*nS-1) = 0d0
      RPA_matrix(2*nS-1,:) = 0d0
      RPA_matrix(:,2*nS) = 0d0
      RPA_matrix(2*nS,:) = 0d0
    end if

    print  *,"U"
    call complex_matout(2*nS,2*nS,complex_U)
    print  *,"Vt"
    call complex_matout(2*nS,2*nS,complex_Vt)
    RPA_matrix         = matmul(matmul(complex_U,RPA_matrix),complex_Vt)
    print  *,"Approx RPA matrix"
    call complex_matout(2*nS,2*nS,RPA_matrix)
    
    TDA_matrix = cmplx(0d0,0d0,kind=8)
    TDA_matrix(1:nS,1:nS)           =   Aph
    TDA_matrix(nS+1:2*nS,nS+1:2*nS) = - Aph
    
    TDA_matrix = matmul(matmul(transpose(conjg(complex_U)),TDA_matrix)&
                                   ,transpose(conjg(complex_Vt)))
    if(abs(D(2*nS))<thresh_reg) then
      TDA_matrix(:,2*nS-1) = 0d0
      TDA_matrix(2*nS-1,:) = 0d0
      TDA_matrix(:,2*nS) = 0d0
      TDA_matrix(2*nS,:) = 0d0
    end if
    TDA_matrix         = matmul(matmul(complex_U,TDA_matrix),complex_Vt)
    Aph                = TDA_matrix(1:nS,1:nS)
    
    !print *, "SVD Regularisation of RPA problem"
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Singularvalue decomp  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !! SVD of RPA matrix
    !allocate(RPA_matrix_real(2*nS,2*nS),U(2*nS,2*nS),Vt(2*nS,2*nS),D(2*nS),&
    !         TDA_matrix_real(2*nS,2*nS))
    !
    !RPA_matrix_real = real(RPA_matrix)
    !TDA_matrix_real = 0d0
    !TDA_matrix_real(1:nS,1:nS)           =   real(Aph)
    !TDA_matrix_real(nS+1:2*nS,nS+1:2*nS) = - real(Aph)

    !print*,"RPA matrix"
    !call matout(2*nS,2*nS,RPA_matrix_real)
    !print*,"TDA matrix"
    !call matout(2*nS,2*nS,TDA_matrix_real)
    !
    !call svd(2*nS,RPA_matrix_real,U,D,Vt)
    !
    !print*, "Singular values of RPA matrix"
    !call vecout(2*nS,D)
    !print*, "U(:,2nS)"
    !call vecout(2*nS,U(:,2*nS))
    !print*, "Vt(:,2nS)"
    !call vecout(2*nS,Vt(2*nS,:))
    !print*, "U(:,2nS-1)"
    !call vecout(2*nS,U(:,2*nS-1))
    !print*, "Vt(:,2nS-1)"
    !call vecout(2*nS,Vt(2*nS-1,:))
    !
    !if(D(2*nS)<thresh_reg) then
    !  
    !  ! Approx TDA matrix
    !  TDA_matrix_real = matmul(matmul(transpose(U),TDA_matrix_real)&
    !                                 ,transpose(Vt))
    !  TDA_matrix_real(:,2*nS) = 0d0
    !  TDA_matrix_real(2*nS,:) = 0d0
    !  TDA_matrix_real         = matmul(matmul(U,TDA_matrix_real),Vt)
    !  Aph                     = cmplx(TDA_matrix_real(1:nS,1:nS),0d0,kind=8)

    !  

    !  ! Approx RPA matrix
    !  D(2*nS) = 0d0 
    !  call AD(2*nS,U,D)
    !  RPA_matrix_real = matmul(U,Vt)
    !  RPA_matrix = cmplx(RPA_matrix_real,0d0,kind=8)
    !  
    !  

    !  print*,"Approx to RPA matrix"
    !  call matout(2*nS,2*nS,RPA_matrix_real)
    !  print*,"Approx to TDA matrix"
    !  call matout(2*nS,2*nS,TDA_matrix_real)
    !  
    !  deallocate(RPA_matrix_real,TDA_matrix_real,U,Vt,D)
    !  
    !  !!! SVD of TDA problem
    !  allocate(RPA_matrix_real(2*nS,2*nS),U(2*nS,2*nS),Vt(2*nS,2*nS),D(2*nS))
    !  RPA_matrix_real                       =  0d0
    !  RPA_matrix_real(1:nS,1:nS)            =  real(Aph)
    !  RPA_matrix_real(nS+1:2*nS,nS+1:2*nS)  = -real(Aph)
    !  call svd(2*nS,RPA_matrix_real,U,D,Vt)
    !  
    !  print*, "Singular values of Aph matrix"
    !  call vecout(2*nS,D)
    !  print*, "U(:,2nS)"
    !  call vecout(2*nS,U(:,2*nS))
    !  print*, "Vt(:,2nS)"
    !  call vecout(2*nS,Vt(2*nS,:))
    !  print*, "U(:,2nS-1)"
    !  call vecout(2*nS,U(:,2*nS-1))
    !  print*, "Vt(:,2nS-1)"
    !  call vecout(2*nS,Vt(2*nS-1,:))
    !  
    !  deallocate(RPA_matrix_real,U,Vt,D)
    !
    ! end if 


    !!!!!!!!!!!!!!!!!!!!!!!!!!

    !RPA_matrix(1:nS,1:nS)           = Aph(:,:)
    !RPA_matrix(1:nS,nS+1:2*nS)      = Bph(:,:)
    !RPA_matrix(nS+1:2*nS,1:nS)      = -Bph(:,:)
    !RPA_matrix(nS+1:2*nS,nS+1:2*nS) = -Aph(:,:)
    
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
    
    !! Print t amplitudes
    !allocate(X_inv(nS,nS))
    !call complex_inverse_matrix(nS,RPA_matrix(1:nS,1:nS),X_inv)
    !call complex_matout(nS,nS,matmul(RPA_matrix(nS+1:2*nS,1:nS),X_inv))
    !deallocate(X_inv)
    
    deallocate(RPA_matrix,OmOmminus)
  
  end if

  ! Compute the RPA correlation energy

  EcRPA = (0.5d0,0.d0)*(sum(Om) - complex_trace_matrix(nS,Aph))
  print *, "1/2 sum Om ", 0.5d0*sum(Om), "1/2 Tr(Aph)", 0.5d0*complex_trace_matrix(nS,Aph)
end subroutine 
