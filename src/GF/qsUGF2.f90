subroutine qsUGF2(maxSCF,thresh,max_diis,BSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,regularize, &
                 nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,EUHF,S,X,T,V,Hc,ERI_AO,               & 
                 ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_AO,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF)

! Perform an unrestricted quasiparticle self-consistent GF2 calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: BSE
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
  double precision              :: Ec(nsp)
  double precision              :: EqsGF2
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: Conv
  double precision              :: rcond(nspin)
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: F_diis(:,:,:)
  double precision,allocatable  :: c(:,:,:)
  double precision,allocatable  :: cp(:,:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: eGF2(:,:)
  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: K(:,:,:)
  double precision,allocatable  :: SigC(:,:,:)
  double precision,allocatable  :: SigCp(:,:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: error(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************************'
  write(*,*)'| Self-consistent unrestricted qsGF2 calculation |'
  write(*,*)'**************************************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsUGF2 !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(eGF2(nBas,nspin),eOld(nBas,nspin),c(nBas,nBas,nspin),cp(nBas,nBas,nspin),P(nBas,nBas,nspin),F(nBas,nBas,nspin), &
           Fp(nBas,nBas,nspin),J(nBas,nBas,nspin),K(nBas,nBas,nspin),SigC(nBas,nBas,nspin),SigCp(nBas,nBas,nspin),         &
           Z(nBas,nspin),error(nBas,nBas,nspin),error_diis(nBasSq,max_diis,nspin),                  &
           F_diis(nBasSq,max_diis,nspin))

! Initialization
  
  nSCF              = -1
  n_diis            = 0
  ispin             = 1
  Conv              = 1d0
  P(:,:,:)          = PHF(:,:,:)
  eGF2(:,:)         = eHF(:,:)
  c(:,:,:)          = cHF(:,:,:)
  F_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0
  rcond(:)          = 0d0

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

    ! 4-index transform for (aa|aa) block

    call AOtoMO_integral_transform(1,1,1,1,nBas,c,ERI_AO,ERI_aaaa)

    ! 4-index transform for (aa|bb) block

    call AOtoMO_integral_transform(1,1,2,2,nBas,c,ERI_AO,ERI_aabb)

    ! 4-index transform for (bb|bb) block

    call AOtoMO_integral_transform(2,2,2,2,nBas,c,ERI_AO,ERI_bbbb)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(regularize) then 

      call UGF2_reg_self_energy(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGF2,SigC,Z)

    else

      call UGF2_self_energy(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGF2,SigC,Z)

    end if

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    do is=1,nspin
      SigC(:,:,is) = 0.5d0*(SigC(:,:,is) + transpose(SigC(:,:,is)))
    end do

    do is=1,nspin
      call MOtoAO_transform(nBas,S,c(:,:,is),SigC(:,:,is),SigCp(:,:,is))
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
      call diagonalize_matrix(nBas,cp(:,:,is),eGF2(:,is))
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

    Conv = maxval(abs(eGF2(:,:) - eOld(:,:)))
    eOld(:,:) = eGF2(:,:)

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
      Ex(is) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,is),K(:,:,is)))
    end do

    ! Correlation energy

    call UMP2(nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EqsGF2,eGF2,Ec)

    ! Total energy

    EqsGF2 = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:)) + sum(Ec(:))

    !------------------------------------------------------------------------
    ! Print results
    !------------------------------------------------------------------------

    call dipole_moment(nBas,P(:,:,1)+P(:,:,2),nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsUGF2(nBas,nO,nSCF,Conv,thresh,eHF,eGF2,c,P,S,T,V,J,K,ENuc,ET,EV,EJ,Ex,Ec,EqsGF2,SigCp,Z,dipole)

  enddo
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

  endif

! Deallocate memory

  deallocate(cp,P,F,Fp,J,K,SigC,SigCp,Z,error,error_diis,F_diis)

! Perform BSE calculation

  if(BSE) then

    print*,'!!! BSE2 NYI for qsUGF2 !!!'

  end if

end subroutine 
