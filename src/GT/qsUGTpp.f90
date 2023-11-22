subroutine qsUGTpp(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE, &
                   TDA_T,TDA,dBSE,dTDA,spin_conserved,spin_flip,&
                   eta,regularize,nBas,nC,nO,nV,nR,nS,nNuc,ZNuc,rNuc,ENuc,EUHF,S,X,T,V,Hc,ERI_AO,ERI_aaaa,&
                   ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GT calculation

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
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip 
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart) 

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin) 
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
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
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart) 

! Local variables
  integer                       :: nSCF
  integer                       :: nBasSq
  double precision              :: dipole(ncart)
  integer                       :: n_diis
  double precision              :: rcond(nspin)
  double precision,external     :: trace_matrix
  double precision              :: Conv
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: Ex(nspin)
  double precision              :: EqsGT
  integer                       :: ispin,is
  integer                       :: iblock 
  integer                       :: nH_sc,nH_sf,nHaa,nHab,nHbb
  integer                       :: nP_sc,nP_sf,nPaa,nPab,nPbb
  double precision              :: EcRPA(nspin),Ecaa,Ecbb
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM(nspin)
  double precision,allocatable  :: Om1ab(:),Om1aa(:),Om1bb(:)
  double precision,allocatable  :: X1ab(:,:),X1aa(:,:),X1bb(:,:)
  double precision,allocatable  :: Y1ab(:,:),Y1aa(:,:),Y1bb(:,:)
  double precision,allocatable  :: rho1ab(:,:,:),rho1aa(:,:,:),rho1bb(:,:,:)
  double precision,allocatable  :: Om2ab(:),Om2aa(:),Om2bb(:)
  double precision,allocatable  :: X2ab(:,:),X2aa(:,:),X2bb(:,:)
  double precision,allocatable  :: Y2ab(:,:),Y2aa(:,:),Y2bb(:,:)
  double precision,allocatable  :: rho2ab(:,:,:),rho2aa(:,:,:),rho2bb(:,:,:)
  double precision,allocatable  :: c(:,:,:)
  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: K(:,:,:) 
  double precision,allocatable  :: SigT(:,:,:)
  double precision,allocatable  :: SigTp(:,:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: eGT(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: e_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)
  double precision,allocatable  :: error(:,:,:)
  

! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Unrestricted evGTpp Calculation *'
  write(*,*)'***********************************'
  write(*,*)

! Dimensions of the pp-URPA linear reponse matrices

  nPaa = nV(1)*(nV(1)-1)/2
  nPbb = nV(2)*(nV(2)-1)/2

  nHaa = nO(1)*(nO(1)-1)/2;
  nHbb = nO(2)*(nO(2)-1)/2;

  nPab = nV(1)*nV(2)
  nHab = nO(1)*nO(2)

  nP_sc = nPab
  nH_sc = nHab

  nP_sf = nPaa + nPbb
  nH_sf = nHaa + nHbb

  nBasSq = nBas*nBas
  

! Memory allocation

  allocate(SigT(nBas,nbas,nspin),SigTp(nBas,nbas,nspin), &
           Z(nBas,nspin),eGT(nBas,nspin),eOld(nBas,nspin), &
           error_diis(nBas,max_diis,nspin),e_diis(nBasSq,max_diis,nspin), &
           F_diis(nBasSq,max_diis,nspin),error(nBas,nBas,nspin),&
           c(nBas,nBas,nspin),cp(nBas,nBas,nspin),P(nBas,nBas,nspin),F(nBas,nBas,nspin), &
           Fp(nBas,nBas,nspin),J(nBas,nBas,nspin),K(nBas,nBas,nspin))

  allocate(Om1ab(nPab),X1ab(nPab,nPab),Y1ab(nHab,nPab), & 
           Om2ab(nHab),X2ab(nPab,nHab),Y2ab(nHab,nHab), & 
           rho1ab(nBas,nBas,nPab),rho2ab(nBas,nBas,nHab), & 
           Om1aa(nPaa),X1aa(nPaa,nPaa),Y1aa(nHaa,nPaa), & 
           Om2aa(nHaa),X2aa(nPaa,nHaa),Y2aa(nHaa,nHaa), & 
           rho1aa(nBas,nBas,nPaa),rho2aa(nBas,nBas,nHaa), & 
           Om1bb(nPbb),X1bb(nPbb,nPbb),Y1bb(nHbb,nPbb), &
           Om2bb(nPbb),X2bb(nPbb,nPbb),Y2bb(nHbb,nPbb), &
           rho1bb(nBas,nBas,nPbb),rho2bb(nBas,nBas,nHbb)) 

!Initialization
  
  nSCF              = -1 
  n_diis            = 0
  Conv              = 1d0
  P(:,:,:)          = PHF(:,:,:) 
  e_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0  
  eGT(:,:)          = eHF(:,:)
  eOld(:,:)         = eHF(:,:)
  c(:,:,:)          = cHF(:,:,:)
  Z(:,:)            = 1d0
  rcond(:)          = 0d0 

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

! Increment

    nSCF = nSCF + 1

! Buid Hartree matrix
    do ispin=1,nspin
      call Hartree_matrix_AO_basis(nBas,P(:,:,ispin),ERI_AO(:,:,:,:), &
                                   J(:,:,ispin))
    end do

! Compute exchange part of the self-energy 
    do ispin=1,nspin
      call exchange_matrix_AO_basis(nBas,P(:,:,ispin),ERI_AO(:,:,:,:), &
                                    K(:,:,ispin))
    end do

! AO to MO transformation of two-electron integrals

     ! 4-index transform for (aa|aa) block

    call AOtoMO_ERI_UHF(1,1,nBas,c,ERI_AO,ERI_aaaa)

    ! 4-index transform for (aa|bb) block

    call AOtoMO_ERI_UHF(1,2,nBas,c,ERI_AO,ERI_aabb)

    ! 4-index transform for (bb|bb) block

    call AOtoMO_ERI_UHF(2,2,nBas,c,ERI_AO,ERI_bbbb) 

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

    ispin  = 1
    iblock = 3
! iblock = 1

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPab,nHaa,nHab,nHbb,nHab,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin)) 
  
! EcRPA(ispin) = 1d0*EcRPA(ispin) 

!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

    ispin  = 2
    iblock = 4

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPaa,nHaa,nHab,nHbb,nHaa,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))  
  
    Ecaa = EcRPA(2)

!----------------------------------------------
! beta-beta block
!----------------------------------------------

    ispin  = 2
    iblock = 7

! Compute linear response

    call ppULR(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb,nPbb,nHaa,nHab,nHbb,nHbb,1d0,eGT,ERI_aaaa, &
               ERI_aabb,ERI_bbbb,Om1bb,X1bb,Y1bb,Om2bb,X2bb,Y2bb,EcRPA(ispin))

    Ecbb = EcRPA(2)
    EcRPA(2) = Ecaa + Ecbb
    EcRPA(1) = EcRPA(1) - EcRPA(2)
    EcRPA(2) = 3d0*EcRPA(2)   

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

!alpha-beta block
  
    iblock = 3

    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHab,nPab,ERI_aaaa,ERI_aabb,ERI_bbbb,X1ab,Y1ab, &
                                  rho1ab,X2ab,Y2ab,rho2ab)
!alpha-alpha block

    iblock = 4
  
    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHaa,nPaa,ERI_aaaa,ERI_aabb,ERI_bbbb,X1aa,Y1aa, &
                                  rho1aa,X2aa,Y2aa,rho2aa)

!beta-beta block 
  
    iblock = 7

    call UGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nHbb,nPbb,ERI_aaaa,ERI_aabb,ERI_bbbb,X1bb,Y1bb, &
                                  rho1bb,X2bb,Y2bb,rho2bb)

    call UGTpp_self_energy(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,nPab,nPbb,eGT,Om1aa,Om1ab,Om1bb,&
                           rho1aa,rho1ab,rho1bb,Om2aa,Om2ab,Om2bb,rho2aa,rho2ab,rho2bb,EcGM,SigT,Z)

! Make correlation self-energy Hermitian and transform it back to AO basis

    do ispin=1,nspin
      SigT(:,:,ispin) = 0.5d0*(SigT(:,:,ispin) + transpose(SigT(:,:,ispin)))
    end do

    do ispin=1,nspin
      call MOtoAO(nBas,S,c(:,:,ispin),SigT(:,:,ispin),SigTp(:,:,ispin))
    end do

! Solve the quasi-particle equation

    do ispin=1,nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + K(:,:,ispin) &
                   + SigTp(:,:,ispin)
    end do    

! Compute commutator and convergence criteria

    do ispin=1,nspin
      error_diis(:,:,ispin) = matmul(F(:,:,ispin),matmul(P(:,:,ispin),S(:,:))) &
                         - matmul(matmul(S(:,:),P(:,:,ispin)),F(:,:,ispin))
    end do 

! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    if(minval(rcond(:)) > 1d-7) then
      do ispin=1,nspin
        if(nO(ispin) > 1) call DIIS_extrapolation(rcond(ispin),nBasSq,nBasSq,n_diis, &
                                                  error_diis(:,1:n_diis,ispin), &
                                                  F_diis(:,1:n_diis,ispin),&
                                                  error_diis(:,:,ispin),F(:,:,ispin))
      end do
    else
      n_diis = 0
    end if

! Transform Fock matrix in orthogonal basis

    do ispin=1,nspin
      Fp(:,:,ispin) = matmul(transpose(X(:,:)),matmul(F(:,:,ispin),X(:,:)))
    end do

! Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:,:) = Fp(:,:,:)
    do ispin=1,nspin
      call diagonalize_matrix(nBas,cp(:,:,ispin),eGT(:,ispin))
    end do

! Back-transform eigenvectors in non-orthogonal basis

    do ispin=1,nspin
      c(:,:,ispin) = matmul(X(:,:),cp(:,:,ispin))
    end do

! Back-transform self-energy

    do ispin=1,nspin
      SigTp(:,:,ispin) = matmul(transpose(c(:,:,ispin)),matmul(SigTp(:,:,ispin),c(:,:,ispin)))
    end do

! Compute density matrix 

    do ispin=1,nspin
      P(:,:,ispin) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
    end do

! Save quasiparticles energy for next cycle

    Conv = maxval(abs(eGT(:,:) - eOld(:,:)))
    eOld(:,:) = eGT(:,:)

!------------------------------------------------------------------------
!   Compute total energy
!------------------------------------------------------------------------

    ! Kinetic energy

    do ispin=1,nspin
      ET(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),T(:,:)))
    end do

    ! Potential energy

    do ispin=1,nspin
      EV(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),V(:,:)))
    end do

! Hartree energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
    EJ(2) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2))) &
          + 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,1)))
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

    ! Exchange energy

    do ispin=1,nspin
      Ex(ispin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin),K(:,:,ispin)))
    end do
write(*,*) 'EcGM', EcGM(1)
    ! Total energy

    EqsGT = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:))

! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsUGT(nBas,nO,nSCF,Conv,thresh,eHF,eGT,c,SigTp,Z,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGT,dipole)

  enddo
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

  endif

! Free memory

  deallocate(Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab, &
             Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, &
             Om1bb,X1bb,Y1bb,Om2bb,X2bb,Y2bb,rho1bb,rho2bb)

  deallocate(c,cp,P,F,Fp,J,K,SigT,SigTp,Z,error,error_diis,F_diis) 

! Testing zone
  
  if(dotest) then
  
    call dump_test_value('U','qsGTpp correlation energy',sum(EcRPA))
    call dump_test_value('U','qsGTpp HOMOa energy',eGT(nO(1),1))
    call dump_test_value('U','qsGTpp LUMOa energy',eGT(nO(1)+1,1))
    call dump_test_value('U','qsGTpp HOMOa energy',eGT(nO(2),2))
    call dump_test_value('U','qsGTpp LUMOa energy',eGT(nO(2)+1,2))
  
  end if

end subroutine 
