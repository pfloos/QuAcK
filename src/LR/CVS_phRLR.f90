subroutine CVS_phRLR(TDA,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

! Compute linear response for unrestricted formalism

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: Aph(nSt,nSt)
  double precision,intent(in)   :: Bph(nSt,nSt)
  
! Local variables

  double precision,external     :: trace_matrix
  double precision,allocatable  :: ApB(:,:)
  double precision,allocatable  :: Aphtilde(:,:)
  double precision,allocatable  :: AmB(:,:)
  double precision,allocatable  :: XpYTDA(:,:)
  double precision,allocatable  :: OmTDA(:)
  complex*16,allocatable        :: complex_Om(:)
  double precision              :: eF
  integer                       :: i

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: Om(nSt)
  double precision,intent(out)  :: XpY(nSt,nSt)
  double precision,intent(out)  :: XmY(nSt,nSt)

! Memory allocation

  allocate(ApB(nSt,nSt),AmB(nSt,nSt),Aphtilde(nSt,nSt))

! Tamm-Dancoff approximation


  if(TDA) then

    XpY(:,:) = Aph(:,:)
    call diagonalize_matrix(nSt,XpY,Om)
    XpY(:,:) = transpose(XpY(:,:))
    XmY(:,:) = XpY(:,:)

    ! Compute the RPA correlation energy
    EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))
  
  else
     
    ! Build (\tilde{A} +B ) (\tilde{A} - B ) - 4*eF*\tilde{A} 
    eF  = 0d0
    Om  = 1d0
    Aphtilde = Aph
    call add_diagonal_matrix(nSt,2*eF*Om,Aphtilde)

    ApB(:,:) = Aphtilde(:,:) + Bph(:,:) 
    AmB(:,:) = Aphtilde(:,:) - Bph(:,:)

    ApB      = matmul(ApB,AmB)
    ApB(:,:) = ApB(:,:) - 4*eF*Aphtilde(:,:)
    call add_diagonal_matrix(nSt,4*eF**2*Om,ApB)

    ! Diagonalize linear response matrix
    call diagonalize_general_matrix(nSt,ApB,Om,XmY)
    
    ! Deal with complex eigenvalue and then stop program
    if(any(Om < 0d0)) then
      print*,"Found complex eigenvalue in Linear Response ! Use complex code !"
      allocate(complex_Om(nSt))
      complex_Om = cmplx(Om,0d0,kind=8)
      complex_Om = sqrt(complex_Om)
      print*,"Excitation energies:"
      call complex_vecout(nSt,complex_Om)
      EcRPA      = 0.5d0*(sum(real(complex_Om)) - trace_matrix(nSt,Aph))
      print*, " Re(EcRPA)  = ", EcRPA 
      print*, "|Im(EcRPA)| = ", 0.5d0*sum(aimag(complex_Om))
      !print *, "Discard imaginary eigenvalue"
      !Om    = real(complex_Om)
      !Om(1) = 0d0
      !
      !allocate(OmTDA(nSt),XpYTDA(nSt,nSt))
    
      !! Compute TDA
      !OmTDA = 0d0
      !XpYTDA(:,:) = Aph(:,:)
      !call diagonalize_matrix(nSt,XpYTDA,OmTDA)
      !print *, "Discarded TDA eigenvalue:", OmTDA(1)
      !OmTDA(1) = 0d0
      !EcRPA = 0.5d0 * (sum(Om) - sum(OmTDA))
      !deallocate(OmTDA,XpYTDA,complex_Om)
      !print *, "Attention: Transition vectors are not computed !"
      stop
    
    else
      
      Om = sqrt(Om)
      
      ! Get XpY via (X+Y) =(\tilde{A}-B - 2*eF*Id)(X-Y)Omega^-1
      XpY = matmul(AmB,XmY) - 2*eF*XmY
      call AD(nSt,XpY,1/Om)
       
      ! Compute Overlap and assign ex-/deexcitations
      ApB = matmul(transpose(XmY),XpY)
      do i=1,nSt
        if(ApB(i,i)<0d0) then
          Om(i)    = - Om(i)
          ! Correct column for right sign of Om 
          XpY(:,i) = - XpY(:,i)
        end if
      enddo
      
      
      call sort_eigenvalues_vec_vec(nSt,Om,XmY,XpY)
      
      ! Orthonormalize
      ApB = matmul(transpose(XmY),XpY)
      call orthogonalize_matrix(1,nSt,ApB,AmB)
      XmY = matmul(XmY,AmB)
      XpY = matmul(XpY,AmB)
      
      ! Transpose XmY and XpY because Quack wants this
      XmY = transpose(XmY) 
      XpY = transpose(XpY)
      
      ! Compute the RPA correlation energy

      EcRPA = 0.5d0*(sum(Om) - trace_matrix(nSt,Aph))

    end if  
  end if

  
  deallocate(ApB,AmB,Aphtilde)

end subroutine 
