subroutine qsUGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
                eta,doSRG,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,ERI_AO,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
                dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)

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
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

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
  double precision,intent(inout):: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA
  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: ixyz
  integer                       :: is
  integer                       :: n_diis
  integer                       :: nSa,nSb,nSt
  double precision              :: dipole(ncart)

  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: EK(nspin)
  double precision              :: EcRPA(nspin)
  double precision              :: EcGM(nspin)
  double precision              :: EqsGW
  double precision              :: EcBSE(nspin)
  double precision              :: Conv
  double precision              :: rcond(nspin)
  double precision,external     :: trace_matrix
  double precision,allocatable  :: err_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:,:)
  double precision,allocatable  :: c(:,:,:)
  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: eGW(:,:)
  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: K(:,:,:)
  double precision,allocatable  :: SigC(:,:,:)
  double precision,allocatable  :: SigCp(:,:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: err(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Unrestricted qsGW Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsUGW !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas
  dRPA = .true.

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! SRG regularization

  if(doSRG) then

    write(*,*) '*** SRG regularized qsGW scheme ***'
    write(*,*)

  end if

! Memory allocation

  nSa = nS(1)
  nSb = nS(2)
  nSt = nSa + nSb

  allocate(Aph(nSt,nSt),Bph(nSt,nSt),eGW(nBas,nspin),c(nBas,nBas,nspin),cp(nBas,nBas,nspin),P(nBas,nBas,nspin), &
           F(nBas,nBas,nspin),Fp(nBas,nBas,nspin),J(nBas,nBas,nspin),K(nBas,nBas,nspin),                        & 
           SigC(nBas,nBas,nspin),SigCp(nBas,nBas,nspin),Z(nBas,nspin),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt),        &  
           rho(nBas,nBas,nSt,nspin),err(nBas,nBas,nspin),err_diis(nBasSq,max_diis,nspin),F_diis(nBasSq,max_diis,nspin))

! Initialization
  
  nSCF            = -1
  n_diis          = 0
  ispin           = 1
  Conv            = 1d0
  P(:,:,:)        = PHF(:,:,:)
  eGW(:,:)        = eHF(:,:)
  c(:,:,:)        = cHF(:,:,:)
  F_diis(:,:,:)   = 0d0
  err_diis(:,:,:) = 0d0
  rcond(:)        = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Hartree matrix

    do is=1,nspin
      call Hartree_matrix_AO_basis(nBas,P(:,:,is),ERI_AO(:,:,:,:),J(:,:,is))
    end do

    ! Compute exchange part of the self-energy 

    do is=1,nspin
      call exchange_matrix_AO_basis(nBas,P(:,:,is),ERI_AO(:,:,:,:),K(:,:,is))
    end do

    !--------------------------------------------------
    ! AO to MO transformation of two-electron integrals
    !--------------------------------------------------
 
    do ixyz=1,ncart
        call AOtoMO(nBas,nBas,c(:,:,1),dipole_int_AO(:,:,ixyz),dipole_int_aa(:,:,ixyz))
        call AOtoMO(nBas,nBas,c(:,:,2),dipole_int_AO(:,:,ixyz),dipole_int_bb(:,:,ixyz))
    end do

    ! 4-index transform for (aa|aa) block

    call AOtoMO_ERI_UHF(1,1,nBas,c,ERI_AO,ERI_aaaa)

    ! 4-index transform for (aa|bb) block

    call AOtoMO_ERI_UHF(1,2,nBas,c,ERI_AO,ERI_aabb)

    ! 4-index transform for (bb|bb) block

    call AOtoMO_ERI_UHF(2,2,nBas,c,ERI_AO,ERI_bbbb)

    ! Compute linear response

    call phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
    
    call phULR(TDA_W,nSa,nSb,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)

    !----------------------!
    ! Excitation densities !
    !----------------------!

    call UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(doSRG) then
      call UGW_SRG_self_energy(nBas,nC,nO,nV,nR,nSt,eGW,Om,rho,SigC,Z,EcGM)
    else
      call UGW_self_energy(eta,nBas,nC,nO,nV,nR,nSt,eGW,Om,rho,SigC,Z,EcGM)
    end if


    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    do is=1,nspin
      SigC(:,:,is) = 0.5d0*(SigC(:,:,is) + transpose(SigC(:,:,is)))
    end do

    do is=1,nspin
      call MOtoAO(nBas,nBas,S,c(:,:,is),SigC(:,:,is),SigCp(:,:,is))
    end do
 
    ! Solve the quasi-particle equation

    do is=1,nspin
      F(:,:,is) = Hc(:,:) + J(:,:,is) + J(:,:,mod(is,2)+1) + K(:,:,is) + SigCp(:,:,is)
    end do

   ! Check convergence 

    do is=1,nspin
      err(:,:,is) = matmul(F(:,:,is),matmul(P(:,:,is),S(:,:))) - matmul(matmul(S(:,:),P(:,:,is)),F(:,:,is))
    end do

    if(nSCF > 1) Conv = maxval(abs(err))

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      do is=1,nspin
        if(nO(is) > 1) call DIIS_extrapolation(rcond(is),nBasSq,nBasSq,n_diis,err_diis(:,1:n_diis,is), &
                                               F_diis(:,1:n_diis,is),err(:,:,is),F(:,:,is))
      end do

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
      call AOtoMO(nBas,nBas,c(:,:,is),SigCp(:,:,is),SigC(:,:,is))
    end do

    ! Compute density matrix 

    do is=1,nspin
      P(:,:,is) = matmul(c(:,1:nO(is),is),transpose(c(:,1:nO(is),is)))
    end do

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

    ! Hartree energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
    EJ(2) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2))) &
          + 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,1)))
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

    ! Exchange energy

    do is=1,nspin
      EK(is) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,is),K(:,:,is)))
    end do

    ! Total energy

    EqsGW = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(EK(:))

    !------------------------------------------------------------------------
    ! Print results
    !------------------------------------------------------------------------

    call dipole_moment(nBas,P(:,:,1)+P(:,:,2),nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsUGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,S,ENuc,ET,EV,EJ,EK,EcGM,EcRPA(ispin),EqsGW,SigCp,Z,dipole)

  end do
!------------------------------------------------------------------------
! End main loop
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

! Deallocate memory

  deallocate(cp,P,F,Fp,J,K,SigC,SigCp,Z,Om,XpY,XmY,rho,err,err_diis,F_diis)

! Perform BSE calculation

  if(BSE) then

    call UGW_phBSE(TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                   S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,c,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 0.5d0*EcBSE(2)

    else

      EcBSE(2) = 0.0d0

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@UHF correlation energy (spin-conserved) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@UHF correlation energy (spin-flip)      = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@UHF correlation energy                  = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@UHF total       energy                  = ',ENuc + EqsGW + sum(EcBSE),' au'
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

      call UGW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip, &
                       eta,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eGW,eGW,EcRPA)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@UHF correlation energy (spin-conserved) = ',EcRPA(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@UHF correlation energy (spin-flip)      = ',EcRPA(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@UHF correlation energy                  = ',sum(EcRPA),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@UHF total       energy                  = ',ENuc + EqsGW + sum(EcRPA),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

! Testing zone
  
  if(dotest) then
  
    call dump_test_value('U','qsGW correlation energy',EcRPA)
    call dump_test_value('U','qsGW HOMOa energy',eGW(nO(1),1))
    call dump_test_value('U','qsGW LUMOa energy',eGW(nO(1)+1,1))
    call dump_test_value('U','qsGW HOMOa energy',eGW(nO(2),2))
    call dump_test_value('U','qsGW LUMOa energy',eGW(nO(2)+1,2))

  end if

end subroutine 
