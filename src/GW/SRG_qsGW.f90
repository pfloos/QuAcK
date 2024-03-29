subroutine SRG_qsGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,BSE2,TDA_W,TDA,          & 
                    dBSE,dTDA,singlet,triplet,eta,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,ERHF, &
                    S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: BSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_MO(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: ixyz
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: Ex
  double precision              :: EqsGW
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  double precision              :: tao,tao1,tao2,tsrg,tsrg1,tsrg2,tlr,tlr1,tlr2,t1,t2,tt,tmo1,tmo2,tmo,tex,tex1,tex2
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)

  logical                       :: dRPA = .true.
  logical                       :: print_W = .true.
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigCp(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error(:,:)

  double precision,parameter    :: flow = 500d0

! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Restricted SRG-qsGW Calculation *'
  write(*,*)'***********************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGW !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

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

  allocate(eGW(nBas),eOld(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas),  &
           J(nBas,nBas),K(nBas,nBas),SigC(nBas,nBas),SigCp(nBas,nBas),Z(nBas),Aph(nS,nS),Bph(nS,nS), & 
           Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS),error(nBas,nBas),error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Initialization
  
  nSCF            = -1
  n_diis          = 0
  ispin           = 1
  Conv            = 1d0
  P(:,:)          = PHF(:,:)
  eGW(:)          = eHF(:)
  eOld(:)         = eHF(:)
  c(:,:)          = cHF(:,:)
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  rcond           = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Hartree matrix
    call wall_time(t1)
    call Hartree_matrix_AO_basis(nBas,P,ERI_AO,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)
    call wall_time(t2)
    tt=tt+t2-t1

    ! AO to MO transformation of two-electron integrals

    call wall_time(tao1)
    
    do ixyz=1,ncart             
      call AOtoMO(nBas,cHF,dipole_int_AO(:,:,ixyz),dipole_int_MO(:,:,ixyz))
    end do  

    call AOtoMO_ERI_RHF(nBas,c,ERI_AO,ERI_MO)

    call wall_time(tao2)

    tao = tao + tao2 -tao1

    ! Compute linear response

    call wall_time(tlr1)

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
    if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

    call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    call wall_time(tlr2)

    tlr = tlr + tlr2 -tlr1

    if(print_W) call print_excitation_energies('phRPA@RGW','singlet',nS,Om)

    ! Compute correlation part of the self-energy 

    call wall_time(tex1)
    
    call GW_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY,rho)

    call wall_time(tex2)
    tex=tex+tex2-tex1

    call wall_time(tsrg1)
    call SRG_self_energy(flow,nBas,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)

    call wall_time(tsrg2)

    tsrg = tsrg + tsrg2 -tsrg1

    ! Make correlation self-energy Hermitian and transform it back to AO basis

    call wall_time(tmo1)
    call MOtoAO(nBas,S,c,SigC,SigCp)
    call wall_time(tmo2)
    tmo = tmo + tmo2 - tmo1
   ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)

    ! Compute commutator and convergence criteria

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    ! DIIS extrapolation 

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,eGW)
    c = matmul(X,cp)

    call AOtoMO(nBas,c,SigCp,SigC)

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Save quasiparticles energy for next cycle

    Conv = maxval(abs(error))
    eOld(:) = eGW(:)

    !------------------------------------------------------------------------
    !   Compute total energy
    !------------------------------------------------------------------------

    ! Kinetic energy

    ET = trace_matrix(nBas,matmul(P,T))

    ! Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    Ex = 0.25d0*trace_matrix(nBas,matmul(P,K))

    ! Total energy

    EqsGW = ET + EV + EJ + Ex 

    ! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsRGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,SigC,Z,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGW,dipole)

  end do
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

 end if

 print *, "Wall time for Fock and exchange build", tt
 print *, "Wall Time for AO to MO", tao
 print *, "Wall Time for LR", tlr
 print *, "Wall Time for excitation density", tex
 print *, "Wall Time for SRG", tsrg
 print *, "Wall time MO to AO Sigma", tmo

! Cumulant expansion

  call RGWC(dotest,eta,nBas,nC,nO,nV,nR,nS,Om,rho,eHF,eGW,eGW,Z)

! Deallocate memory

  deallocate(c,cp,P,F,Fp,J,K,SigC,Z,Om,XpY,XmY,rho,error,error_diis,F_diis)

! Perform BSE calculation

  if(BSE) then

    call GW_phBSE(BSE2,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(2)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call GW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (singlet) =',EcBSE(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (triplet) =',EcBSE(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy           =',EcBSE(1) + EcBSE(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine 
