subroutine UKS(x_rung,x_DFA,c_rung,c_DFA,nEns,wEns,nCC,aCC,nGrid,weight,maxSCF,thresh,max_diis,guess_type,mix,level_shift, &
               nNuc,ZNuc,rNuc,ENuc,nBas,AO,dAO,S,T,V,Hc,ERI,dipole_int,X,occnum,Cx_choice,doNcentered,Ew,eKS,c,Pw,Vxc)

! Perform unrestricted Kohn-Sham calculation for ensembles

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x_rung,c_rung
  integer,intent(in)            :: x_DFA,c_DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  logical,intent(in)            :: mix
  double precision,intent(in)   :: level_shift
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: occnum(nBas,nspin,nEns)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered

! Local variables

  integer                       :: xc_rung
  logical                       :: LDA_centered = .false.
  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: n_diis
  integer                       :: nO(nspin,nEns)
  double precision              :: Conv
  double precision              :: rcond(nspin)
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EH(nsp)
  double precision              :: Ex(nspin)
  double precision              :: Ec(nsp)
  double precision              :: dipole(ncart)

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

  double precision,allocatable  :: rhow(:,:)
  double precision,allocatable  :: drhow(:,:,:)
  double precision              :: nEl(nspin)

  double precision,allocatable  :: P(:,:,:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: drho(:,:,:,:)

  integer                       :: ispin,iEns,iBas

! Output variables

  double precision,intent(out)  :: Ew
  double precision,intent(out)  :: eKS(nBas,nspin)
  double precision,intent(out)  :: Pw(nBas,nBas,nspin)
  double precision,intent(out)  :: c(nBas,nBas,nspin)
  double precision,intent(out)  :: Vxc(nBas,nspin)

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

  xc_rung = max(x_rung,c_rung)

! Memory allocation

  allocate(cp(nBas,nBas,nspin),J(nBas,nBas,nspin),F(nBas,nBas,nspin),Fp(nBas,nBas,nspin),      & 
           Fx(nBas,nBas,nspin),FxHF(nBas,nBas,nspin),Fc(nBas,nBas,nspin),err(nBas,nBas,nspin), &
           rhow(nGrid,nspin),drhow(ncart,nGrid,nspin),                                         &
           err_diis(nBasSq,max_diis,nspin),F_diis(nBasSq,max_diis,nspin),                      &
           P(nBas,nBas,nspin,nEns),rho(nGrid,nspin,nEns),drho(ncart,nGrid,nspin,nEns))

! Guess coefficients and eigenvalues

  nO(:,:) = 0
  do iEns=1,nEns
    do ispin=1,nspin
      nO(ispin,iEns) = int(sum(occnum(:,ispin,iEns)))
    end do
  end do

  do ispin=1,nspin
    call mo_guess(nBas,guess_type,S,Hc,X,c(:,:,ispin))
  end do

! Mix guess for UHF solution in singlet states

  if(mix) then
    write(*,*) '!! guess mixing disabled in UKS !!'
    write(*,*) 
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
  rcond(:)        = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','E(KS)','|','Ex(KS)','|','Ec(KS)','|','Conv','|','nEl','|'
  write(*,*)'------------------------------------------------------------------------------------------'
  
  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!------------------------------------------------------------------------
!   Compute density matrix 
!------------------------------------------------------------------------

    call density_matrix(nBas,nEns,c(:,:,:),P(:,:,:,:),occnum(:,:,:))
 
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

    if(xc_rung > 1) then 

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

!------------------------------------------------------------------------
!   Compute Hxc potential and Fock operator
!------------------------------------------------------------------------

!   Compute Hartree potential 

    do ispin=1,nspin
      call hartree_potential(nBas,Pw(:,:,ispin),ERI,J(:,:,ispin))
    end do

!   Compute exchange potential

    do ispin=1,nspin
      call exchange_potential(x_rung,x_DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas, &
                                           Pw(:,:,ispin),ERI,AO,dAO,rhow(:,ispin),drhow(:,:,ispin),       &
                                           Cx_choice,doNcentered,Fx(:,:,ispin),FxHF(:,:,ispin))
    end do

!   Compute correlation potential

    call correlation_potential(c_rung,c_DFA,nEns,wEns,nGrid,weight,nBas,AO,dAO,rhow,drhow,Fc)

!   Build Fock operator

    do ispin=1,nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + Fx(:,:,ispin) + Fc(:,:,ispin)
    end do

!   Check convergence 

    do ispin=1,nspin
      err(:,:,ispin) = matmul(F(:,:,ispin),matmul(Pw(:,:,ispin),S(:,:))) - matmul(matmul(S(:,:),Pw(:,:,ispin)),F(:,:,ispin))
    end do

    if(nSCF > 1) Conv = maxval(abs(err(:,:,:)))
    
!   DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    do ispin=1,nspin

      if(rcond(ispin) > 1d-15) then

        call DIIS_extrapolation(rcond(ispin),nBasSq,nBasSq,n_diis, & 
                                err_diis(:,:,ispin),F_diis(:,:,ispin),err(:,:,ispin),F(:,:,ispin))
      else 

        n_diis = 0

      end if

    end do

!   Level-shifting

    if(level_shift > 0d0 .and. Conv > thresh) then 

      do ispin=1,nspin
        call level_shifting(level_shift,nBas,maxval(nO(ispin,:)),S,c,F(:,:,ispin))
      end do

    end if

!  Transform Fock matrix in orthogonal basis

    do ispin=1,nspin
      Fp(:,:,ispin) = matmul(transpose(X(:,:)),matmul(F(:,:,ispin),X(:,:)))
    end do

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:,:) = Fp(:,:,:)
    do ispin=1,nspin
      call diagonalize_matrix(nBas,cp(:,:,ispin),eKS(:,ispin))
    end do
    
!   Back-transform eigenvectors in non-orthogonal basis

    do ispin=1,nspin
      c(:,:,ispin) = matmul(X(:,:),cp(:,:,ispin))
    end do

!------------------------------------------------------------------------
!   Compute KS energy
!------------------------------------------------------------------------

!   Kinetic energy

    do ispin=1,nspin
      ET(ispin) = trace_matrix(nBas,matmul(Pw(:,:,ispin),T(:,:)))
    end do

!   Potential energy

    do ispin=1,nspin
      EV(ispin) = trace_matrix(nBas,matmul(Pw(:,:,ispin),V(:,:)))
    end do

!   Hartree energy

    call hartree_energy(nBas,Pw,J,EH)

!   Exchange energy

    do ispin=1,nspin
      call exchange_energy(x_rung,x_DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas, &
                                        Pw(:,:,ispin),FxHF(:,:,ispin),rhow(:,ispin),drhow(:,:,ispin),  &
                                        Cx_choice,doNcentered,Ex(ispin))
    end do

!   Correlation energy

    call correlation_energy(c_rung,c_DFA,nEns,wEns,nGrid,weight,rhow,drhow,Ec)

!   Total energy

    Ew = sum(ET(:)) + sum(EV(:)) + sum(EH(:)) + sum(Ex(:)) + sum(Ec(:))

!   Check the grid accuracy by computing the number of electrons 

    do ispin=1,nspin
      nEl(ispin) = electron_number(nGrid,weight,rhow(:,ispin))
    end do

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',Ew + ENuc,'|',sum(Ex(:)),'|',sum(Ec(:)),'|',Conv,'|',sum(nEl(:)),'|'

  end do
  write(*,*)'------------------------------------------------------------------------------------------'

!  print*,'Ensemble energy:',Ew + ENuc,'au'
  

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

  call dipole_moment(nBas,Pw(:,:,1)+Pw(:,:,2),nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_UKS(nBas,nEns,occnum,S,wEns,eKS,c,ENuc,ET,EV,EH,Ex,Ec,Ew,dipole)

! Compute Vxc for post-HF calculations

  call xc_potential(nBas,c,Fx,Fc,Vxc)

!------------------------------------------------------------------------
! Compute individual energies from ensemble energy
!------------------------------------------------------------------------

  call individual_energy(x_rung,x_DFA,c_rung,c_DFA,LDA_centered,nEns,wEns,nCC,aCC,nGrid,weight,nBas,      &
                                      AO,dAO,T,V,ERI,ENuc,eKS,Pw,rhow,drhow,J,Fx,FxHF,Fc,P,rho,drho,occnum,Cx_choice,doNcentered,Ew)

end subroutine UKS