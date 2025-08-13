subroutine G_qsParquet(TDAeh,TDApp,max_diis_1b,max_diis_2b,eta_1b,eta_2b,ENuc,max_it_1b,conv_1b,max_it_2b,conv_2b, & 
                    nOrb,nOrb2,nC,nO,nV,nR,nS,EGHF,PHF,cHF,eHF,Ov,Or,T,V,Hc,ERI_AO,ERI_MO)

! Parquet approximation with quasiparticle self-consistency based on spin orbitals

  implicit none
  include 'parameters.h'

! Hard-coded parameters

  logical                       :: print_phLR = .false. ! Print the eh two-body excitations at each two-body iterations
  logical                       :: print_ppLR = .false. ! Print the pp two-body excitations at each two-body iterations
  
! Input variables

  logical,intent(in)            :: TDAeh      
  logical,intent(in)            :: TDApp      
  integer,intent(in)            :: max_diis_1b
  integer,intent(in)            :: max_diis_2b
  double precision,intent(in)   :: eta_1b,eta_2b        
  double precision,intent(in)   :: ENuc
  integer,intent(in)            :: max_it_1b,max_it_2b
  double precision,intent(in)   :: conv_1b,conv_2b
  integer,intent(in)            :: nOrb,nOrb2,nC,nO,nV,nR,nS
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: PHF(nOrb2,nOrb2)
  double precision,intent(in)   :: cHF(nOrb2,nOrb2)
  double precision,intent(in)   :: eHF(nOrb2)
  double precision,intent(in)   :: Ov(nOrb,nOrb)
  double precision,intent(in)   :: Or(nOrb,nOrb)
  double precision,intent(in)   :: T(nOrb,nOrb)
  double precision,intent(in)   :: V(nOrb,nOrb)
  double precision,intent(in)   :: Hc(nOrb,nOrb)
  double precision,intent(in)   :: ERI_AO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(inout):: ERI_MO(nOrb2,nOrb2,nOrb2,nOrb2)

! Local variables

  integer                       :: n_it_1b,n_it_2b
  double precision              :: err_1b
  double precision              :: err_2b
  double precision              :: err_eh, err_pp
  double precision              :: err_eig_eh,err_eig_pp,err_eig_hh,err_eig_ee
  double precision              :: start_t,end_t,tt
  double precision              :: start_1b,end_1b,t_1b
  double precision              :: start_2b,end_2b,t_2b

  integer                       :: nOO,nVV,nOrb2_Sq

  ! eh BSE
  double precision              :: Ec_eh
  double precision,allocatable  :: Aph(:,:), Bph(:,:)
  double precision,allocatable  :: XpY(:,:), XmY(:,:)
  double precision,allocatable  :: eh_Om(:), old_eh_Om(:)
  double precision,allocatable  :: eh_Gam_A(:,:),eh_Gam_B(:,:)
  ! pp BSE
  double precision              :: Ec_pp
  double precision,allocatable  :: Bpp(:,:), Cpp(:,:), Dpp(:,:)
  double precision,allocatable  :: X1(:,:),Y1(:,:)
  double precision,allocatable  :: ee_Om(:), old_ee_Om(:)
  double precision,allocatable  :: X2(:,:),Y2(:,:)
  double precision,allocatable  :: hh_Om(:), old_hh_Om(:)
  double precision,allocatable  :: pp_Gam_B(:,:),pp_Gam_C(:,:),pp_Gam_D(:,:)
  ! Effective integrals
  double precision,allocatable  :: eh_rho(:,:,:), ee_rho(:,:,:), hh_rho(:,:,:)
  ! Reducible kernels
  double precision,allocatable  :: eh_Phi(:,:,:,:), pp_Phi(:,:,:,:)
  double precision,allocatable  :: old_eh_Phi(:,:,:,:), old_pp_Phi(:,:,:,:)
  ! One-body quantities
  double precision,allocatable  :: Ca(:,:),Cb(:,:)
  double precision,allocatable  :: ERI_tmp(:,:,:,:)
  double precision,allocatable  :: Jaa(:,:),Jbb(:,:)
  double precision,allocatable  :: Kaa(:,:),Kab(:,:),Kba(:,:),Kbb(:,:)
  double precision,allocatable  :: Faa(:,:),Fab(:,:),Fba(:,:),Fbb(:,:)
  double precision,allocatable  :: Paa(:,:),Pab(:,:),Pba(:,:),Pbb(:,:)
  
  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: Cp(:,:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: eQP(:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigC_AO(:,:)
  double precision,allocatable  :: Z(:)
  double precision              :: EcGM
  ! DIIS
  integer                       :: n_diis_1b,n_diis_2b
  double precision              :: rcond_1b,rcond_2b
  double precision,allocatable  :: err_diis_1b(:,:),err_diis_2b(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: err_F(:,:)
  double precision,allocatable  :: Phi_diis(:,:)
  double precision,allocatable  :: err(:)
  double precision,allocatable  :: Phi(:)
  double precision              :: alpha

  integer                       :: pp,q,r,ss,pqrs

  double precision              :: mem = 0d0
  double precision              :: dp_in_GB = 8d0/(1024d0**3)

! Output variables
! None
  
! Useful parameters
  nOO = nO*(nO - 1)/2
  nVV = nV*(nV - 1)/2
  nOrb2_Sq = nOrb2 * nOrb2

! Start
 
  write(*,*)
  write(*,*)'*************************************'
  write(*,*)'* Generalized qsParquet Calculation *'
  write(*,*)'*************************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsParquet !!'
  write(*,*)

! Print parameters

  write(*,*)'---------------------------------------------------------------'
  write(*,*)' Parquet parameters for one-body and two-body self-consistency '
  write(*,*)'---------------------------------------------------------------'
  write(*,'(1X,A50,1X,I5)')    'Maximum number of one-body iteration:',max_it_1b
  write(*,'(1X,A50,1X,E10.5)') 'Convergence threshold for one-body energies:',conv_1b
  write(*,'(1X,A50,1X,E10.5)') 'Strength of SRG one-body regularization:',eta_1b
  write(*,'(1X,A50,1X,E10.5)') 'Strength of SRG two-body regularization:',eta_2b
  write(*,'(1X,A50,1X,I5)')    'Maximum length of DIIS expansion:',max_diis_1b
  write(*,*)'---------------------------------------------------------------'
  write(*,'(1X,A50,1X,I5)')    'Maximum number of two-body iteration:',max_it_2b
  write(*,'(1X,A50,1X,E10.5)') 'Convergence threshold for two-body energies:',conv_2b
  write(*,'(1X,A50,1X,L5)')    'TDA for eh excitation energies?',TDAeh
  write(*,'(1X,A50,1X,L5)')    'TDA for pp excitation energies?',TDApp
  write(*,'(1X,A50,1X,I5)')    'Maximum length of DIIS expansion:',max_diis_2b
  write(*,*)'---------------------------------------------------------------'
  write(*,*)
    
! DIIS for two-body part

  allocate(err_diis_2b(2*nOrb2**4,max_diis_2b),Phi_diis(2*nOrb2**4,max_diis_2b))
  allocate(err(2*nOrb2**4),Phi(2*nOrb2**4))

  mem = mem + size(err_diis_2b) + size(Phi_diis) + size(err) + size(Phi)
  write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

! DIIS for one-body part        

  allocate(err_diis_1b(nOrb2_Sq,max_diis_1b),F_diis(nOrb2_Sq,max_diis_1b),err_F(nOrb2,nOrb2))

  mem = mem + size(err_diis_1b) + size(F_diis) + size(err_F)
  write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet = ',mem*dp_in_GB,' GB'

  rcond_1b  = 0d0
  n_diis_1b = 0
  err_diis_1b(:,:) = 0d0
  F_diis(:,:)    = 0d0
  
! Memory allocation for Parquet quantity
  
  allocate(eQP(nOrb2))

  mem = mem + size(eQP) 

  allocate(old_eh_Om(nS),old_ee_Om(nVV),old_hh_Om(nOO))
  allocate(eh_rho(nOrb2,nOrb2,nS),ee_rho(nOrb2,nOrb2,nVV),hh_rho(nOrb2,nOrb2,nOO))
  allocate(old_eh_Phi(nOrb2,nOrb2,nOrb2,nOrb2),old_pp_Phi(nOrb2,nOrb2,nOrb2,nOrb2))

  mem = mem + size(old_eh_Om) + size(old_ee_Om) + size(old_hh_Om)
  mem = mem + size(eh_rho) + size(ee_rho) + size(hh_rho)
  mem = mem + size(old_eh_Phi) + size(old_pp_Phi)

! Memory allocation for qs self-consistency

  allocate(Ca(nOrb,nOrb2),Cb(nOrb,nOrb2),Jaa(nOrb,nOrb),Jbb(nOrb,nOrb),   &
           Kaa(nOrb,nOrb),Kab(nOrb,nOrb),Kba(nOrb,nOrb),Kbb(nOrb,nOrb),   &
           Faa(nOrb,nOrb),Fab(nOrb,nOrb),Fba(nOrb,nOrb),Fbb(nOrb,nOrb),   &
           Paa(nOrb,nOrb),Pab(nOrb,nOrb),Pba(nOrb,nOrb),Pbb(nOrb,nOrb),   &
           P(nOrb2,nOrb2),H(nOrb2,nOrb2),S(nOrb2,nOrb2),X(nOrb2,nOrb2),   &
           C(nOrb2,nOrb2),Cp(nOrb2,nOrb2),F(nOrb2,nOrb2),Fp(nOrb2,nOrb2))

  mem = mem + size(Ca)  + size(Cb)  + size(Jaa) + size(Jbb)
  mem = mem + size(Kaa) + size(Kab) + size(Kba) + size(Kbb)
  mem = mem + size(Faa) + size(Fab) + size(Fba) + size(Fbb)
  mem = mem + size(Paa) + size(Pab) + size(Pba) + size(Pbb)
  mem = mem + size(P)   + size(H)   + size(S)   + size(X)
  mem = mem + size(C)   + size(Cp)  + size(F)   + size(Fp)
  write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

! Initialization

  n_it_1b = 0
  err_1b  = 1d0

  EcGM = 0d0
  Ec_eh = 0d0 ! TODO: This should be removed once the GM is coded.
  Ec_pp = 0d0
  
  eQP(:)  = eHF(:)  
  P(:,:)        = PHF(:,:)
  C(:,:)        = cHF(:,:)

! Construct super overlap matrix

  S(      :     ,       :    ) = 0d0
  S(     1:nOrb ,     1:nOrb ) = Ov(1:nOrb,1:nOrb)
  S(nOrb+1:nOrb2,nOrb+1:nOrb2) = Ov(1:nOrb,1:nOrb)

! Construct super orthogonalization matrix

  X(      :     ,       :    ) = 0d0
  X(     1:nOrb ,     1:nOrb ) = Or(1:nOrb,1:nOrb)
  X(nOrb+1:nOrb2,nOrb+1:nOrb2) = Or(1:nOrb,1:nOrb)

! Construct super orthogonalization matrix

  H(      :     ,       :    ) = 0d0
  H(     1:nOrb ,     1:nOrb ) = Hc(1:nOrb,1:nOrb)
  H(nOrb+1:nOrb2,nOrb+1:nOrb2) = Hc(1:nOrb,1:nOrb)

! Construct super density matrix

  P(:,:) = matmul(C(:,1:nO),transpose(C(:,1:nO)))

  Paa(:,:) = P(     1:nOrb ,     1:nOrb )
  Pab(:,:) = P(     1:nOrb ,nOrb+1:nOrb2)
  Pba(:,:) = P(nOrb+1:nOrb2,     1:nOrb )
  Pbb(:,:) = P(nOrb+1:nOrb2,nOrb+1:nOrb2)    
  
  !-----------------------------------------!
  ! Main loop for one-body self-consistency !
  !-----------------------------------------!

  do while(err_1b > conv_1b .and. n_it_1b < max_it_1b)

    n_it_1b = n_it_1b + 1
    call wall_time(start_1b)

    write(*,*)
    write(*,*)'====================================='
    write(*,'(1X,A30,1X,I4)') 'One-body iteration #',n_it_1b
    write(*,*)'====================================='
    write(*,*)

!   Buid Hartree matrix

    call Hartree_matrix_AO_basis(nOrb,Paa,ERI_AO,Jaa)
    call Hartree_matrix_AO_basis(nOrb,Pbb,ERI_AO,Jbb)

!   Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nOrb,Paa,ERI_AO,Kaa)
    call exchange_matrix_AO_basis(nOrb,Pba,ERI_AO,Kab)
    call exchange_matrix_AO_basis(nOrb,Pab,ERI_AO,Kba)
    call exchange_matrix_AO_basis(nOrb,Pbb,ERI_AO,Kbb)

!   Build individual Fock matrices

    Faa(:,:) = Hc(:,:) + Jaa(:,:) + Jbb(:,:) + Kaa(:,:) 
    Fab(:,:) =                               + Kab(:,:)
    Fba(:,:) =                               + Kba(:,:)
    Fbb(:,:) = Hc(:,:) + Jbb(:,:) + Jaa(:,:) + Kbb(:,:)

!   Build super Fock matrix

    F(     1:nOrb ,     1:nOrb ) = Faa(1:nOrb,1:nOrb)
    F(     1:nOrb ,nOrb+1:nOrb2) = Fab(1:nOrb,1:nOrb)
    F(nOrb+1:nOrb2,     1:nOrb ) = Fba(1:nOrb,1:nOrb)
    F(nOrb+1:nOrb2,nOrb+1:nOrb2) = Fbb(1:nOrb,1:nOrb)

!   AO to MO transformation of two-electron integrals

    allocate(ERI_tmp(nOrb2,nOrb2,nOrb2,nOrb2))

    mem = mem + size(ERI_tmp)
    write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'
    
    Ca(:,:) = C(1:nOrb,1:nOrb2)
    Cb(:,:) = C(nOrb+1:nOrb2,1:nOrb2)

!TODO I copy pasted this from qsGGW but that would be great to avoid a third copy of the ERI
    call AOtoMO_ERI_GHF(nOrb,nOrb2,Ca,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nOrb,nOrb2,Ca,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nOrb,nOrb2,Cb,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nOrb,nOrb2,Cb,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    deallocate(ERI_tmp)
    mem = mem - size(ERI_tmp)
    write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'
    
!   Initialization
    
    n_it_2b = 0 
    err_2b  = 1d0
    ! We cannot use the quantity from the previous one-body iteration because the orbitals changed
    eh_rho(:,:,:) = 0d0
    ee_rho(:,:,:) = 0d0
    hh_rho(:,:,:) = 0d0

    old_eh_Om(:) = 0d0
    old_ee_Om(:) = 0d0
    old_hh_Om(:) = 0d0
    
    old_eh_Phi(:,:,:,:) = 0d0
    old_pp_Phi(:,:,:,:) = 0d0
    ! Same for DIIS
    rcond_2b  = 1d0
    n_diis_2b = 0
    err_diis_2b(:,:) = 0d0
    Phi_diis(:,:)    = 0d0

    !-----------------------------------------!
    ! Main loop for two-body self-consistency !
    !-----------------------------------------!
    
    do while(err_2b > conv_2b .and. n_it_2b < max_it_2b)

      n_it_2b = n_it_2b + 1
      call wall_time(start_2b)

      write(*,*)' ***********************************'
      write(*,'(1X,A30,1X,I4)') 'Two-body iteration #',n_it_2b
      write(*,*)' ***********************************'
      write(*,*)

      !-----------------!
      !    eh channel   !
      !-----------------!

      write(*,*) 'Diagonalizing eh BSE problem...'

      allocate(Aph(nS,nS),Bph(nS,nS),eh_Om(nS),XpY(nS,nS),XmY(nS,nS),eh_Gam_A(nS,nS),eh_Gam_B(nS,nS))

      mem = mem + size(Aph) + size(Bph) + size(eh_Om) + size(XpY) + size(XmY) + size(eh_Gam_A) + size(eh_Gam_B)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      Aph(:,:) = 0d0
      Bph(:,:) = 0d0

      call wall_time(start_t)

                     call phGLR_A(.false.,nOrb2,nC,nO,nV,nR,nS,1d0,eQP,ERI_MO,Aph)
      if(.not.TDAeh) call phGLR_B(.false.,nOrb2,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)
      
      if(n_it_1b == 1 .and. n_it_2b == 1) then

        eh_Gam_A(:,:) = 0d0
        eh_Gam_B(:,:) = 0d0

      else

                       call G_eh_Gamma_A(nOrb2,nC,nO,nR,nS,old_eh_Phi,old_pp_Phi,eh_Gam_A)
        if(.not.TDAeh) call G_eh_Gamma_B(nOrb2,nC,nO,nR,nS,old_eh_Phi,old_pp_Phi,eh_Gam_B)
        
      end if
      
      Aph(:,:) = Aph(:,:) + eh_Gam_A(:,:)
      Bph(:,:) = Bph(:,:) + eh_Gam_B(:,:) 
      
      call phGLR(TDAeh,nS,Aph,Bph,Ec_eh,eh_Om,XpY,XmY)

 !     call matout(nS,nS,XpY)

      call wall_time(end_t)

      tt = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for phBSE problem =',tt,' seconds'
      write(*,*)

      if(print_phLR) call print_excitation_energies('phBSE@Parquet','eh generalized',nS,eh_Om)

      err_eig_eh = maxval(abs(old_eh_Om - eh_Om))

      deallocate(Aph,Bph,eh_Gam_A,eh_Gam_B)

      mem = mem - size(Aph) - size(Bph) - size(eh_Gam_A) - size(eh_Gam_B)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      !-----------------!
      !    pp channel   !
      !-----------------!

      write(*,*) 'Diagonalizing pp BSE problem...'

      allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO),   &
               ee_Om(nVV),X1(nVV,nVV),Y1(nOO,nVV), & 
               hh_Om(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
               pp_Gam_B(nVV,nOO),pp_Gam_C(nVV,nVV),pp_Gam_D(nOO,nOO))

      mem = mem + size(Bpp) + size(Cpp) + size(Dpp) &
                + size(ee_Om) + size(X1) + size(Y1) &
                + size(hh_Om) + size(X2) + size(Y2) &
                + size(pp_Gam_B) + size(pp_Gam_C) + size(pp_Gam_D)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'


      Bpp(:,:) = 0d0
      Cpp(:,:) = 0d0
      Dpp(:,:) = 0d0

      call wall_time(start_t)
      if(.not.TDApp) call ppGLR_B(nOrb2,nC,nO,nV,nR,nOO,nVV,1d0,ERI_MO,Bpp)
                     call ppGLR_C(nOrb2,nC,nO,nV,nR,nVV,1d0,eQP,ERI_MO,Cpp)
                     call ppGLR_D(nOrb2,nC,nO,nV,nR,nOO,1d0,eQP,ERI_MO,Dpp)

      if(n_it_1b == 1 .and. n_it_2b == 1) then

        pp_Gam_B(:,:) = 0d0
        pp_Gam_C(:,:) = 0d0
        pp_Gam_D(:,:) = 0d0

      else

        if(.not.TDApp) call G_pp_Gamma_B(nOrb2,nC,nO,nR,nOO,nVV,old_eh_Phi,pp_Gam_B)
                       call G_pp_Gamma_C(nOrb2,nO,nR,nVV,old_eh_Phi,pp_Gam_C)
                       call G_pp_Gamma_D(nOrb2,nC,nO,nOO,old_eh_Phi,pp_Gam_D)

      end if
                   
      Bpp(:,:) = Bpp(:,:) + pp_Gam_B(:,:)
      Cpp(:,:) = Cpp(:,:) + pp_Gam_C(:,:)
      Dpp(:,:) = Dpp(:,:) + pp_Gam_D(:,:)
      
      call ppGLR(TDApp,nOO,nVV,Bpp,Cpp,Dpp,ee_Om,X1,Y1,hh_Om,X2,Y2,Ec_pp)
      call wall_time(end_t)
      tt = end_t - start_t

      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for ppBSE problem =',tt,' seconds'
      write(*,*)

      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p generalized',nVV,ee_Om)
      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h generalized',nOO,hh_Om)

      err_eig_ee = maxval(abs(old_ee_Om - ee_Om))
      err_eig_hh = maxval(abs(old_hh_Om - hh_Om))
      err_eig_pp = max(err_eig_ee,err_eig_hh)

      deallocate(Bpp,Cpp,Dpp,pp_Gam_B,pp_Gam_C,pp_Gam_D)

      mem = mem - size(Bpp) - size(Cpp) - size(Dpp) &
                - size(pp_Gam_B) - size(pp_Gam_C) - size(pp_Gam_D)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      !----------!
      ! Updating !
      !----------!

      old_eh_Om(:) = eh_Om(:)
      old_ee_Om(:) = ee_Om(:)
      old_hh_Om(:) = hh_Om(:)

      deallocate(eh_Om,ee_Om,hh_Om)

      mem = mem - size(eh_Om) - size(ee_Om) - size(hh_Om)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'
      
      !----------------------------!
      ! Compute screened integrals !
      !----------------------------!

      ! Free memory
      deallocate(eh_rho,ee_rho,hh_rho)

      mem = mem - size(eh_rho) - size(ee_rho) - size(hh_rho)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      ! TODO Once we will compute the blocks of kernel starting from the 4-tensors we can move the freeing up
      ! Memory allocation
      allocate(eh_rho(nOrb2,nOrb2,nS))
      allocate(ee_rho(nOrb2,nOrb2,nVV),hh_rho(nOrb2,nOrb2,nOO))

      mem = mem + size(eh_rho) + size(ee_rho) + size(hh_rho)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      ! Build singlet eh integrals
      write(*,*) 'Computing eh screened integrals...'

      call wall_time(start_t)
      call G_eh_screened_integral(nOrb2,nC,nO,nR,nS,ERI_MO,old_eh_Phi,old_pp_Phi,XpY,XmY,eh_rho)
      call wall_time(end_t)
      tt = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for eh integrals =',tt,' seconds'
      write(*,*)

      ! Done with eigenvectors and kernel

      deallocate(XpY,XmY)

      mem = mem - size(XpY) - size(XmY)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      ! Build singlet pp integrals
      write(*,*) 'Computing pp screened integrals...'

      call wall_time(start_t)
      call G_pp_screened_integral(nOrb2,nC,nO,nR,nOO,nVV,ERI_MO,old_eh_Phi,X1,Y1,ee_rho,X2,Y2,hh_rho)
      call wall_time(end_t)
      tt = end_t - start_t
      ! Done with eigenvectors and kernel

      deallocate(X1,Y1,X2,Y2)

      mem = mem - size(X1) - size(Y1) - size(X2) - size(Y2)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for pp integrals =',tt,' seconds'
      write(*,*)

      !----------------------------!
      ! Compute reducible kernels  !
      !----------------------------!

      ! Memory allocation
      allocate(eh_Phi(nOrb2,nOrb2,nOrb2,nOrb2))
      allocate(pp_Phi(nOrb2,nOrb2,nOrb2,nOrb2))

      mem = mem + size(eh_Phi) + size(pp_Phi)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

      ! Build eh reducible kernels
      write(*,*) 'Computing eh reducible kernel...'

      call wall_time(start_t)
      call G_eh_Phi(eta_2b,nOrb2,nC,nR,nS,old_eh_Om,eh_rho,eh_Phi)
      call wall_time(end_t)
      tt = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for eh reducible kernel =',tt,' seconds'
      write(*,*)

      ! Build pp reducible kernels
      write(*,*) 'Computing pp reducible kernel...'

      call wall_time(start_t)
      call G_pp_Phi(eta_2b,nOrb2,nC,nR,nOO,nVV,old_ee_Om,ee_rho,old_hh_Om,hh_rho,pp_Phi)
      call wall_time(end_t)
      tt = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for pp reducible kernel =',tt,' seconds'
      write(*,*)

      err_eh = maxval(abs(eh_Phi - old_eh_Phi))
      err_pp = maxval(abs(pp_Phi - old_pp_Phi))

      ! alpha = 0.5d0
      ! eh_Phi(:,:,:,:) = alpha * eh_Phi(:,:,:,:) + (1d0 - alpha) * old_eh_Phi(:,:,:,:)
      ! pp_Phi(:,:,:,:) = alpha * pp_Phi(:,:,:,:) + (1d0 - alpha) * old_pp_Phi(:,:,:,:)

      !--------------------!
      ! DIIS extrapolation !
      !--------------------!

      pqrs = 0
      do pp=1,nOrb2
        do q=1,nOrb2
          do r=1,nOrb2
            do ss=1,nOrb2
              pqrs = pqrs + 1

              err(        pqrs) = eh_Phi(pp,q,r,ss) - old_eh_Phi(pp,q,r,ss)
              err(nOrb2**4+pqrs) = pp_Phi(pp,q,r,ss) - old_pp_Phi(pp,q,r,ss)

              Phi(        pqrs) = eh_Phi(pp,q,r,ss)
              Phi(nOrb2**4+pqrs) = pp_Phi(pp,q,r,ss)

            end do
          end do
        end do
      end do

      if(max_diis_2b > 1) then
     
        n_diis_2b = min(n_diis_2b+1,max_diis_2b)
        call DIIS_extrapolation(rcond_2b,2*nOrb2**4,2*nOrb2**4,n_diis_2b,err_diis_2b,Phi_diis,err,Phi)
     
      end if

      pqrs = 0
      do pp=1,nOrb2
        do q=1,nOrb2
          do r=1,nOrb2
            do ss=1,nOrb2
              pqrs = pqrs + 1

              eh_Phi(pp,q,r,ss) = Phi(        pqrs)
              pp_Phi(pp,q,r,ss) = Phi(nOrb2**4+pqrs) 

            end do
          end do
        end do
      end do

      old_eh_Phi(:,:,:,:) = eh_Phi(:,:,:,:)
      old_pp_Phi(:,:,:,:) = pp_Phi(:,:,:,:)
      
      ! Free memory

      deallocate(eh_Phi,pp_Phi)

      mem = mem - size(eh_Phi) - size(pp_Phi)
      write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'
      
      write(*,*) '------------------------------------------------'
      write(*,*) '    Two-body (frequency/kernel) convergence     '
      write(*,*) '------------------------------------------------'
      write(*,'(1X,A24,F10.6,1X,A1,1X,F10.6)')'Error for eh channel = ',err_eig_eh,'/',err_eh
      write(*,'(1X,A24,F10.6,1X,A1,1X,F10.6)')'Error for pp channel = ',err_eig_pp,'/',err_pp
      write(*,*) '------------------------------------------------'
      write(*,*)
      
      ! Convergence criteria
      err_2b = max(err_eh,err_pp)
      
      call wall_time(end_2b)
      t_2b = end_2b - start_2b
      write(*,'(1X,A44,1X,I4,A2,F9.3,A8)') 'Wall time for two-body iteration #',n_it_2b,' =',t_2b,' seconds'
      write(*,*)

    end do
    !---------------------------------------------!
    ! End main loop for two-body self-consistency !
    !---------------------------------------------!

    ! Did it actually converge?

    if(n_it_2b == max_it_2b) then

      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'             Two-body convergence failed            '
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      !stop

    else

      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'             Two-body convergence success           '
      write(*,*)'****************************************************'
      write(*,*)

      call print_excitation_energies('phBSE@Parquet','1h1p',nS,old_eh_Om)
      call print_excitation_energies('ppBSE@Parquet','2p',nVV,old_ee_Om)
      call print_excitation_energies('ppBSE@Parquet','2h',nOO,old_hh_Om)

    end if

    allocate(Z(nOrb2),SigC(nOrb2,nOrb2),SigC_AO(nOrb2,nOrb2)) 
    
    mem = mem + size(Z) + size(SigC) + size(SigC_AO)
    write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

    write(*,*) 'Computing self-energy...'
    write(*,*) 
    
    call wall_time(start_t)
    call G_Parquet_self_energy(eta_1b,nOrb2,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI_MO, &
                               eh_rho,old_eh_Om,ee_rho,old_ee_Om,hh_rho,old_hh_Om,SigC,Z)
    call wall_time(end_t)
    tt = end_t - start_t
    write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for self energy =',tt,' seconds'
    write(*,*)

    SigC_AO(:,:) = 0d0
    call MOtoAO(nOrb2,nOrb2,S,C,SigC,SigC_AO)

!   ... and add self-energy
    
    F(:,:) = F(:,:) + SigC_AO(:,:)

!   Compute commutator and convergence criteria

    err_F  = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    err_1b = maxval(abs(err_F))

!   DIIS for one-body part

    if(max_diis_1b > 1) then 
  
      n_diis_1b = min(n_diis_1b+1,max_diis_1b)
      call DIIS_extrapolation(rcond_1b,nOrb2_Sq,nOrb2_Sq,n_diis_1b,err_diis_1b,F_diis,err_F,F)
  
    end if

!  Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X),matmul(F,X))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    Cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nOrb2,Cp,eQP)

!   Back-transform eigenvectors in non-orthogonal basis

    C(:,:) = matmul(X,Cp)

    !call AOtoMO_GHF(nOrb,nOrb2,Ca,Cb,SigC_AO,SigC)
    call AOtoMO(nOrb2,nOrb2,C,SigC_AO,SigC)

!   Form super density matrix

    P(:,:) = matmul(C(:,1:nO),transpose(C(:,1:nO)))

!   Compute individual density matrices

    Paa(:,:) = P(     1:nOrb ,     1:nOrb )
    Pab(:,:) = P(     1:nOrb ,nOrb+1:nOrb2)
    Pba(:,:) = P(nOrb+1:nOrb2,     1:nOrb )
    Pbb(:,:) = P(nOrb+1:nOrb2,nOrb+1:nOrb2)
    
    ! Print for one-body part
    call G_print_qsparquet_1b(nOrb2,nC,nO,nV,nR,eHF,SigC,eQP,Z,n_it_1b,err_1b,ENuc,EGHF,EcGM,Ec_eh,Ec_pp)

    mem = mem - size(Z) - size(SigC)
    write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'

    call wall_time(end_1b)
    t_1b = end_1b - start_1b
    write(*,'(1X,A44,1X,I4,A2,F9.3,A8)') 'Wall time for one-body iteration #',n_it_1b,' =',t_1b,' seconds'

    deallocate(Z,SigC,SigC_AO)
    
    mem = mem - size(Z) - size(SigC) - size(SigC_AO)
    write(*,'(1X,A50,4X,F6.3,A3)') 'Memory usage in GParquet =',mem*dp_in_GB,' GB'
    
  end do 
  !---------------------------------------------!
  ! End main loop for one-body self-consistency !
  !---------------------------------------------!
  
  ! Did it actually converge?
  if(n_it_1b == max_it_1b) then

     write(*,*)
     write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)'             One-body convergence failed            '
     write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)
     ! stop
     
  else
     
     write(*,*)
     write(*,*)'****************************************************'
     write(*,*)'             One-body convergence success           '
     write(*,*)'****************************************************'
     write(*,*)
     
  end if

  ! call G_Parquet_Galitskii_Migdal(eta_1b,nOrb2,nC,nO,nV,nR,nS,nOO,nVV,eOld,ERI_MO, &
  !                              eh_rho,old_eh_Om,ee_rho,old_ee_Om,hh_rho,old_hh_Om,EcGM)
  
end subroutine 
