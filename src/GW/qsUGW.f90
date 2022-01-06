subroutine qsUGW(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,    & 
                G0W,GW0,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO, & 
                nV,nR,nS,EUHF,S,X,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa, &
                dipole_int_bb,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: G0W
  logical,intent(in)            :: GW0
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)

  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: doGWPT = .false.

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: is
  integer                       :: n_diis
  integer                       :: nS_aa,nS_bb,nS_sc
  double precision              :: dipole(ncart)

  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: Ex(nspin)
  double precision              :: EcRPA
  double precision              :: EcGM(nspin)
  double precision              :: EqsGW
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: Conv
  double precision              :: rcond(nspin)
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)
  double precision,allocatable  :: c(:,:,:)
  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: eGW(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: K(:,:,:)
  double precision,allocatable  :: SigC(:,:,:)
  double precision,allocatable  :: SigCp(:,:,:)
  double precision,allocatable  :: SigCm(:,:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: error(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'*************************************************'
  write(*,*)'| Self-consistent unrestricted qsGW calculation |'
  write(*,*)'*************************************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsUGW !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! COHSEX approximation

  if(COHSEX) then 
    write(*,*) 'COHSEX approximation activated!'
    write(*,*)
  end if

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(eGW(nBas,nspin),eOld(nBas,nspin),c(nBas,nBas,nspin),cp(nBas,nBas,nspin),P(nBas,nBas,nspin),F(nBas,nBas,nspin), &
           Fp(nBas,nBas,nspin),J(nBas,nBas,nspin),K(nBas,nBas,nspin),SigC(nBas,nBas,nspin),SigCp(nBas,nBas,nspin),        &
           SigCm(nBas,nBas,nspin),Z(nBas,nspin),OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),                   &
           rho_RPA(nBas,nBas,nS_sc,nspin),error(nBas,nBas,nspin),error_diis(nBasSq,max_diis,nspin),                       & 
           F_diis(nBasSq,max_diis,nspin))

! Initialization
  
  nSCF              = -1
  n_diis            = 0
  ispin             = 1
  Conv              = 1d0
  P(:,:,:)          = PHF(:,:,:)
  eGW(:,:)          = eHF(:,:)
  eOld(:,:)         = eHF(:,:)
  c(:,:,:)          = cHF(:,:,:)
  F_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0
  rcond             = 1d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Coulomb matrix

    do is=1,nspin
      call Coulomb_matrix_AO_basis(nBas,P(:,:,is),ERI_AO(:,:,:,:),J(:,:,is))
    end do

    ! Compute exchange part of the self-energy 

    do is=1,nspin
      call exchange_matrix_AO_basis(nBas,P(:,:,is),ERI_AO(:,:,:,:),K(:,:,is))
    end do

    !--------------------------------------------------
    ! AO to MO transformation of two-electron integrals
    !--------------------------------------------------

    ! 4-index transform for (aa|aa) block

    call AOtoMO_integral_transform(1,1,1,1,nBas,c,ERI_AO,ERI_aaaa)

    ! 4-index transform for (aa|bb) block

    call AOtoMO_integral_transform(1,1,2,2,nBas,c,ERI_AO,ERI_aabb)

    ! 4-index transform for (bb|bb) block

    call AOtoMO_integral_transform(2,2,2,2,nBas,c,ERI_AO,ERI_bbbb)

    ! Compute linear response

    if(.not. GW0 .or. nSCF == 0) then

      call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                        eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

    endif

    !----------------------!
    ! Excitation densities !
    !----------------------!

    call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(G0W) then

      if(regularize) then

        call unrestricted_regularized_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,SigC,EcGM)
        call unrestricted_regularized_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,Z)

      else

        call unrestricted_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,SigC,EcGM)
        call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,Z)

      end if

     else

      if(regularize) then

        call unrestricted_regularized_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,SigC,EcGM)
        call unrestricted_regularized_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,Z)

      else

        call unrestricted_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,SigC,EcGM)
        call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,Z)

      end if

     endif

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    do is=1,nspin
      SigCp(:,:,is) = 0.5d0*(SigC(:,:,is) + transpose(SigC(:,:,is)))
      SigCm(:,:,is) = 0.5d0*(SigC(:,:,is) - transpose(SigC(:,:,is)))
    end do

    do is=1,nspin
      call MOtoAO_transform(nBas,S,c(:,:,is),SigCp(:,:,is))
    end do
 
    ! Solve the quasi-particle equation

    do is=1,nspin
      F(:,:,is) = Hc(:,:) + J(:,:,is) + J(:,:,mod(is,2)+1) + K(:,:,is) + SigCp(:,:,is)
    end do

   ! Check convergence 

    do is=1,nspin
      error(:,:,is) = matmul(F(:,:,is),matmul(P(:,:,is),S(:,:))) - matmul(matmul(S(:,:),P(:,:,is)),F(:,:,is))
    end do

    if(nSCF > 1) conv = maxval(abs(error(:,:,:)))

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    if(minval(rcond(:)) > 1d-7) then
      do is=1,nspin
        if(nO(is) > 1) call DIIS_extrapolation(rcond(is),nBasSq,nBasSq,n_diis,error_diis(:,1:n_diis,is), &
                                                  F_diis(:,1:n_diis,is),error(:,:,is),F(:,:,is))
      end do
    else 
      n_diis = 0
    end if

    ! Transform Fock matrix in orthogonal basis

    do is=1,nspin
      Fp(:,:,is) = matmul(transpose(X(:,:)),matmul(F(:,:,is),X(:,:)))
    end do

    ! Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:,:) = Fp(:,:,:)
    do is=1,nspin
      call diagonalize_matrix(nBas,cp(:,:,is),eGW(:,is))
    end do

    ! Back-transform eigenvectors in non-orthogonal basis

    do is=1,nspin
      c(:,:,is) = matmul(X(:,:),cp(:,:,is))
    end do

    ! Back-transform self-energy

    do is=1,nspin
      SigCp(:,:,is) = matmul(transpose(c(:,:,is)),matmul(SigCp(:,:,is),c(:,:,is)))
    end do

    ! Compute density matrix 

    do is=1,nspin
      P(:,:,is) = matmul(c(:,1:nO(is),is),transpose(c(:,1:nO(is),is)))
    end do

    ! Save quasiparticles energy for next cycle

    Conv = maxval(abs(eGW(:,:) - eOld(:,:)))
    eOld(:,:) = eGW(:,:)

    !------------------------------------------------------------------------
    !   Compute total energy
    !------------------------------------------------------------------------

    ! Kinetic energy

    do is=1,nspin
      ET(is) = trace_matrix(nBas,matmul(P(:,:,is),T(:,:)))
    end do

    ! Potential energy

    do is=1,nspin
      EV(is) = trace_matrix(nBas,matmul(P(:,:,is),V(:,:)))
    end do

    ! Coulomb energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
    EJ(2) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2))) &
          + 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,1)))
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

    ! Exchange energy

    do is=1,nspin
      Ex(is) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,is),K(:,:,is)))
    end do

    ! Total energy

    EqsGW = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:))

    !------------------------------------------------------------------------
    ! Print results
    !------------------------------------------------------------------------

    call dipole_moment(nBas,P(:,:,1)+P(:,:,2),nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsUGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,P,S,T,V,J,K,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGW,SigCp,Z,dipole)

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Compute second-order correction of the Hermitization error

!if(doGWPT) call qsGW_PT(nBas,nC,nO,nV,nR,nS,eGW,SigCm)

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

! Deallocate memory

  deallocate(cp,P,F,Fp,J,K,SigC,SigCp,SigCm,Z,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,error,error_diis,F_diis)

! Perform BSE calculation

  if(BSE) then

    call unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                                     S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,c,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 0.5d0*EcBSE(2)

    else

      EcBSE(2) = 0.0d0

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsUGW correlation energy (spin-conserved) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsUGW correlation energy (spin-flip)      =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsUGW correlation energy                  =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsUGW total energy                        =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '--------------------------------------------------------------'
      write(*,*) ' Adiabatic connection version of BSE@qsUGW correlation energy '
      write(*,*) '--------------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call unrestricted_ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip, &
                              eta,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eGW,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsUGW correlation energy (spin-conserved) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsUGW correlation energy (spin-flip)      =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsUGW correlation energy                  =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsUGW total energy                        =',ENuc + EqsGW + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine qsUGW
