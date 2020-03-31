subroutine GOK_RKS(restart,x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nGrid,weight,maxSCF,thresh, & 
                   max_diis,guess_type,nBas,AO,dAO,nO,nV,S,T,V,Hc,ERI,X,ENuc,Ew,EwGIC,F)

! Perform restricted Kohn-Sham calculation for ensembles

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: restart
  integer,intent(in)            :: x_rung,c_rung
  character(len=12),intent(in)  :: x_DFA,c_DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: maxSCF,max_diis,guess_type
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ENuc

  double precision,intent(inout):: F(nBas,nBas)

! Local variables

  integer                       :: xc_rung
  integer                       :: nSCF,nBasSq
  integer                       :: n_diis
  double precision              :: conv
  double precision              :: rcond
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: Ex
  double precision              :: Ec

  double precision,allocatable  :: eps(:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: Fx(:,:)
  double precision,allocatable  :: FxHF(:,:)
  double precision,allocatable  :: Fc(:,:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,external     :: trace_matrix
  double precision,external     :: electron_number

  double precision,allocatable  :: Pw(:,:)
  double precision,allocatable  :: rhow(:)
  double precision,allocatable  :: drhow(:,:)
  double precision              :: nEl

  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: rho(:,:)
  double precision,allocatable  :: drho(:,:,:)

  double precision              :: E(nEns)
  double precision              :: Om(nEns)

  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: Ew
  double precision,intent(out)  :: EwGIC

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*      Restricted Kohn-Sham calculation        *'
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
  write(*,*) '*                        EXCHANGE RUNG                            *'
  write(*,*) '*******************************************************************'

  call select_rung(x_rung,x_DFA)

! Select rung for correlation

  write(*,*)
  write(*,*) '*******************************************************************'
  write(*,*) '*                       CORRELATION RUNG                          *'
  write(*,*) '*******************************************************************'

  call select_rung(c_rung,c_DFA)

! Overall rung

  xc_rung = max(x_rung,c_rung)

! Memory allocation

  allocate(eps(nBas),c(nBas,nBas),cp(nBas,nBas),              &
           J(nBas,nBas),Fp(nBas,nBas),Fx(nBas,nBas),          & 
           FxHF(nBas,nBas),Fc(nBas,nBas),err(nBas,nBas),      &
           Pw(nBas,nBas),rhow(nGrid),drhow(ncart,nGrid),      &
           err_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis), &
           P(nBas,nBas,nEns),rho(nGrid,nEns),drho(ncart,nGrid,nEns))

! Guess coefficients and eigenvalues

  if(.not. restart) then
    if(guess_type == 1) then  

    cp(:,:) = matmul(transpose(X(:,:)),matmul(Hc(:,:),X(:,:)))
    call diagonalize_matrix(nBas,cp(:,:),eps(:))
    c(:,:) = matmul(X(:,:),cp(:,:))

    else
 
      print*,'Wrong guess option'
      stop
 
    end if

  end if
    
! Initialization

  nSCF = 0
  conv = 1d0

  nEl = 0d0
  Ex  = 0d0
  Ec  = 0d0

  Fx(:,:)   = 0d0
  FxHF(:,:) = 0d0
  Fc(:,:)   = 0d0

  n_diis        = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0

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

    call restricted_density_matrix(nBas,nEns,nO,c(:,:),P(:,:,:))
 
!   Weight-dependent density matrix
    
    Pw(:,:) = 0d0
    do iEns=1,nEns
      Pw(:,:) = Pw(:,:) + wEns(iEns)*P(:,:,iEns)
    end do

!------------------------------------------------------------------------
!   Compute one-electron density and its gradient if necessary
!------------------------------------------------------------------------

    do iEns=1,nEns
      call density(nGrid,nBas,P(:,:,iEns),AO(:,:),rho(:,iEns))
    end do

!   Weight-dependent one-electron density 

    rhow(:) = 0d0
    do iEns=1,nEns
      rhow(:) = rhow(:) + wEns(iEns)*rho(:,iEns) 
    end do

    if(xc_rung > 1) then 

!     Compute gradient of the one-electron density

      do iEns=1,nEns
        call gradient_density(nGrid,nBas,P(:,:,iEns),AO(:,:),dAO(:,:,:),drho(:,:,iEns))
      end do

!     Weight-dependent one-electron density gradient

      drhow(:,:) = 0d0
      do iEns=1,nEns
        drhow(:,:) = drhow(:,:) + wEns(iEns)*drho(:,:,iEns)
      end do

    end if

!   Build Coulomb repulsion

    call hartree_coulomb(nBas,Pw(:,:),ERI(:,:,:,:),J(:,:))

!   Compute exchange potential

    call exchange_potential(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas,Pw(:,:),ERI(:,:,:,:), &
                            AO(:,:),dAO(:,:,:),rhow(:),drhow(:,:),Fx(:,:),FxHF(:,:))

!   Compute correlation potential

    call restricted_correlation_potential(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:), & 
                                          nBas,AO(:,:),dAO(:,:,:),rhow(:),drhow(:,:),Fc(:,:))

!   Build Fock operator

    F(:,:) = Hc(:,:) + J(:,:) + Fx(:,:) + Fc(:,:)

!   Check convergence 

    err(:,:) = matmul(F(:,:),matmul(Pw(:,:),S(:,:))) - matmul(matmul(S(:,:),Pw(:,:)),F(:,:))

    conv = maxval(abs(err(:,:)))
    
!   DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,err_diis(:,:),F_diis(:,:),err(:,:),F(:,:))

!   Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

!   Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X(:,:)),matmul(F(:,:),X(:,:)))

!   Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp(:,:),eps(:))
    
!   Back-transform eigenvectors in non-orthogonal basis

    c(:,:) = matmul(X(:,:),cp(:,:))

!------------------------------------------------------------------------
!   Compute KS energy
!------------------------------------------------------------------------

!  Kinetic energy

    ET = trace_matrix(nBas,matmul(Pw(:,:),T(:,:)))

!  Potential energy

    EV = trace_matrix(nBas,matmul(Pw(:,:),V(:,:)))

!  Coulomb energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(Pw(:,:),J(:,:)))

!   Exchange energy

    call exchange_energy(x_rung,x_DFA,nEns,wEns(:),nGrid,weight(:),nBas, &
                         Pw(:,:),FxHF(:,:),rhow(:),drhow(:,:),Ex)

!   Correlation energy

    call restricted_correlation_energy(c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),Ec)

!   Total energy

    Ew = ET + EV + EJ + Ex + Ec

!   Check the grid accuracy by computing the number of electrons 

    nEl = electron_number(nGrid,weight(:),rhow(:))

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',Ew + ENuc,'|',Ex,'|',Ec,'|',conv,'|',nEl,'|'
 
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

  call print_RKS(nBas,nO,eps(:),c(:,:),ENuc,ET,EV,EJ,Ex,Ec,Ew)

!------------------------------------------------------------------------
! Compute individual energies from ensemble energy
!------------------------------------------------------------------------

  call restricted_individual_energy(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns(:),nGrid,weight(:), &
                                    nBas,nO,nV,T(:,:),V(:,:),ERI(:,:,:,:),ENuc,             & 
                                    eps(:),Pw(:,:),rhow(:),drhow(:,:),J(:,:),P(:,:,:),      & 
                                    rho(:,:),drho(:,:,:),Ew,EwGIC,E(:),Om(:))

end subroutine GOK_RKS
