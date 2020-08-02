subroutine eDFT_UKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,maxSCF,thresh,max_diis,guess_type, &
                   nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,Ew,occnum,Cx_choice)

! Perform unrestricted Kohn-Sham calculation for ensembles

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: maxSCF,max_diis,guess_type
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  integer,intent(in)            :: nO(nspin),nV(nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: occnum(2,2,3)
  integer,intent(in)            :: Cx_choice

! Local variables

  integer                       :: xc_rung
  logical                       :: LDA_centered = .false.
  integer                       :: nSCF,nBasSq
  integer                       :: n_diis
  double precision              :: conv
  double precision              :: rcond(nspin)
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: Ex(nspin)
  double precision              :: Ec(nsp)
  double precision              :: Ew

  double precision,allocatable  :: eps(:,:)
  double precision,allocatable  :: c(:,:,:)
  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: Fx(:,:,:)
  double precision,allocatable  :: FxHF(:,:,:)
  double precision,allocatable  :: Fc(:,:,:)
  double precision,allocatable  :: err(:,:,:)
  double precision,allocatable  :: err_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)
  double precision,external     :: trace_matrix
  double precision,external     :: electron_number

  double precision,allocatable  :: Pw(:,:,:)
  double precision,allocatable  :: rhow(:,:)
  double precision,allocatable  :: drhow(:,:,:)
  double precision              :: nEl(nspin)

  double precision,allocatable  :: P(:,:,:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: drho(:,:,:,:)

  double precision              :: E(nEns)
  double precision              :: Om(nEns)

  integer                       :: ispin,iEns

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*    Unrestricted Kohn-Sham calculation        *'
  write(*,*)'*           *** for ensembles ***              *'
  write(*,*)'************************************************'
  write(*,*)

! Useful stuff

  nBasSq = nBas*nBas

!------------------------------------------------------------------------
! Rung of Jacob's ladder
!------------------------------------------------------------------------

! Select rung for exchange 

  write(*,*)
  write(*,*) '*******************************************************************'
  write(*,*) '*                        Exchange rung                            *'
  write(*,*) '*******************************************************************'

  call select_rung(x_rung,x_DFA)

! Select rung for correlation

  write(*,*)
  write(*,*) '*******************************************************************'
  write(*,*) '*                       Correlation rung                          *'
  write(*,*) '*******************************************************************'

  call select_rung(c_rung,c_DFA)

! Overall rung

  xc_rung = max(x_rung,c_rung)

! Memory allocation

  allocate(eps(nBas,nspin),c(nBas,nBas,nspin),cp(nBas,nBas,nspin),        &
           J(nBas,nBas,nspin),F(nBas,nBas,nspin),Fp(nBas,nBas,nspin),     & 
           Fx(nBas,nBas,nspin),FxHF(nBas,nBas,nspin),Fc(nBas,nBas,nspin),err(nBas,nBas,nspin), &
           Pw(nBas,nBas,nspin),rhow(nGrid,nspin),drhow(ncart,nGrid,nspin), &
           err_diis(nBasSq,max_diis,nspin),F_diis(nBasSq,max_diis,nspin),  &
           P(nBas,nBas,nspin,nEns),rho(nGrid,nspin,nEns),drho(ncart,nGrid,nspin,nEns))

! Guess coefficients and eigenvalues


  if(guess_type == 1) then

    do ispin=1,nspin
      cp(:,:,ispin) = matmul(transpose(X(:,:)),matmul(Hc(:,:),X(:,:)))
      call diagonalize_matrix(nBas,cp(:,:,ispin),eps(:,ispin))
      c(:,:,ispin) = matmul(X(:,:),cp(:,:,ispin))
    end do

  else if(guess_type == 2) then

    do ispin=1,nspin
      call random_number(F(:,:,ispin))
    end do

  else

    print*,'Wrong guess option'
    stop

  end if

! Initialization

  nSCF = 0
  conv = 1d0

  nEl(:) = 0d0
  Ex(:)  = 0d0
  Ec(:)  = 0d0

  Fx(:,:,:)   = 0d0
  FxHF(:,:,:) = 0d0
  Fc(:,:,:)   = 0d0

  n_diis          = 0
  F_diis(:,:,:)   = 0d0
  err_diis(:,:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','E(KS)','|','Ex(KS)','|','Ec(KS)','|','Conv','|','nEl','|'
  write(*,*)'------------------------------------------------------------------------------------------'
  
  do while(conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!------------------------------------------------------------------------
!   Compute density matrix 
!------------------------------------------------------------------------

    call unrestricted_density_matrix(nBas,nEns,nO(:),c(:,:,:),P(:,:,:,:),occnum)
 
!   Weight-dependent density matrix
    
    Pw(:,:,:) = 0d0
    do iEns=1,nEns
      Pw(:,:,:) = Pw(:,:,:) + wEns(iEns)*P(:,:,:,iEns)
    end do

!------------------------------------------------------------------------
!   Compute one-electron density and its gradient if necessary
!------------------------------------------------------------------------

    do ispin=1,nspin
      do iEns=1,nEns
        call density(nGrid,nBas,P(:,:,ispin,iEns),AO(:,:),rho(:,ispin,iEns))
      end do
    end do

!   Weight-dependent one-electron density 

    rhow(:,:) = 0d0
    do iEns=1,nEns
      rhow(:,:) = rhow(:,:) + wEns(iEns)*rho(:,:,iEns) 
    end do

    if(xc_rung > 1 .and. xc_rung /= 666) then 

!     Ground state density 

      do ispin=1,nspin
        do iEns=1,nEns
          call gradient_density(nGrid,nBas,P(:,:,ispin,iEns),AO(:,:),dAO(:,:,:),drho(:,:,ispin,iEns))
        end do
      end do

!     Weight-dependent one-electron density 

      drhow(:,:,:) = 0d0
      do iEns=1,nEns
        drhow(:,:,:) = drhow(:,:,:) + wEns(iEns)*drho(:,:,:,iEns)
      end do

    end if

!   Build Coulomb repulsion

    do ispin=1,nspin
      call hartree_coulomb(nBas,Pw(:,:,ispin),ERI(:,:,:,:),J(:,:,ispin))
    end do

!   Compute exchange potential

    do ispin=1,nspin
      call unrestricted_exchange_potential(x_rung,x_DFA,LDA_centered,nEns,wEns(:),aCC_w1,aCC_w2,nGrid,weight(:),nBas,    &
                                           Pw(:,:,ispin),ERI(:,:,:,:),AO(:,:),dAO(:,:,:),rhow(:,ispin),drhow(:,:,ispin), &
                                           Fx(:,:,ispin),FxHF(:,:,ispin),Cx_choice)
    end do

!   Compute correlation potential

    call unrestricted_correlation_potential(c_rung,c_DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rhow,drhow,Fc)

!   Build Fock operator
    do ispin=1,nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + Fx(:,:,ispin) + Fc(:,:,ispin)
    end do

!   Check convergence 

    do ispin=1,nspin
      err(:,:,ispin) = matmul(F(:,:,ispin),matmul(Pw(:,:,ispin),S(:,:))) - matmul(matmul(S(:,:),Pw(:,:,ispin)),F(:,:,ispin))
    end do

    conv = maxval(abs(err(:,:,:)))
    
!   DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    do ispin=1,nspin
      call DIIS_extrapolation(rcond(ispin),nBasSq,nBasSq,n_diis, & 
                              err_diis(:,:,ispin),F_diis(:,:,ispin),err(:,:,ispin),F(:,:,ispin))
    end do

!   Reset DIIS if required

    if(minval(rcond(:)) < 1d-15) n_diis = 0

!  Transform Fock matrix in orthogonal basis

    do ispin=1,nspin
      Fp(:,:,ispin) = matmul(transpose(X(:,:)),matmul(F(:,:,ispin),X(:,:)))
    end do

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:,:) = Fp(:,:,:)
    do ispin=1,nspin
      call diagonalize_matrix(nBas,cp(:,:,ispin),eps(:,ispin))
    end do
    
!   Back-transform eigenvectors in non-orthogonal basis

    do ispin=1,nspin
      c(:,:,ispin) = matmul(X(:,:),cp(:,:,ispin))
    end do

!------------------------------------------------------------------------
!   Compute KS energy
!------------------------------------------------------------------------

!  Kinetic energy

    do ispin=1,nspin
      ET(ispin) = trace_matrix(nBas,matmul(Pw(:,:,ispin),T(:,:)))
    end do

!  Potential energy

    do ispin=1,nspin
      EV(ispin) = trace_matrix(nBas,matmul(Pw(:,:,ispin),V(:,:)))
    end do

!  Coulomb energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,1)))
    EJ(2) = 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,1),J(:,:,2))) &
          + 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,1))) 
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(Pw(:,:,2),J(:,:,2)))

!   Exchange energy

    do ispin=1,nspin
      call unrestricted_exchange_energy(x_rung,x_DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas, &
                                        Pw(:,:,ispin),FxHF(:,:,ispin),rhow(:,ispin),drhow(:,:,ispin),Ex(ispin),Cx_choice)
    end do

!   Correlation energy

    call unrestricted_correlation_energy(c_rung,c_DFA,nEns,wEns,nGrid,weight,rhow,drhow,Ec)

!   Total energy

    Ew = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:)) + sum(Ec(:))

!   Check the grid accuracy by computing the number of electrons 

    do ispin=1,nspin
      nEl(ispin) = electron_number(nGrid,weight,rhow(:,ispin))
    end do

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',Ew + ENuc,'|',sum(Ex(:)),'|',sum(Ec(:)),'|',conv,'|',sum(nEl(:)),'|'

  end do
  write(*,*)'------------------------------------------------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

! Compute final KS energy

  call print_UKS(nBas,nO,eps,c,ENuc,ET,EV,EJ,Ex,Ec,Ew)

!------------------------------------------------------------------------
! Compute individual energies from ensemble energy
!------------------------------------------------------------------------

  call unrestricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,nBas, &
                                      AO,dAO,nO,nV,T,V,ERI,ENuc,eps,Pw,rhow,drhow,J,Fx,FxHF,Fc,P,rho,drho,Ew,E,Om,occnum,Cx_choice)

end subroutine eDFT_UKS
