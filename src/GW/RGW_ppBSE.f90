subroutine RGW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies at the pp level

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: nOO
  integer                       :: nVV

  integer                       :: p, q, m
  integer                       :: i_data, supp_data_dbl_size, supp_data_int_size
  integer                       :: n_states, n_states_diag

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)

  double precision,allocatable  :: Om1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)

  double precision,allocatable  :: Om2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

  double precision,allocatable  :: KB_sta(:,:)
  double precision,allocatable  :: KC_sta(:,:)
  double precision,allocatable  :: KD_sta(:,:)
  
  integer,         allocatable  :: supp_data_int(:)
  double precision,allocatable  :: supp_data_dbl(:)
  double precision,allocatable  :: Om(:), R(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)


!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS))
 
  call phLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)

  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  deallocate(XpY_RPA,XmY_RPA,Aph,Bph)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    write(*,*) '****************'
    write(*,*) '*** Singlets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 1
    EcBSE(ispin) = 0d0

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV))
    allocate(Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    ! Compute BSE excitation energies

    ! ---
    ! LAPACK
    ! ---

    allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
    allocate(KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

    call RGW_ppBSE_static_kernel_C(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,OmRPA,rho_RPA,KC_sta)
    call RGW_ppBSE_static_kernel_D(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,OmRPA,rho_RPA,KD_sta)
    if(.not.TDA) call RGW_ppBSE_static_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,KB_sta)
    
    call ppLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
    call ppLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
    if(.not.TDA) call ppLR_B(ispin,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))
    deallocate(Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)

    !print*, 'LAPACK:'
    !print*, Om2
    !print*, Om1

    ! ---



    ! ---
    ! Davidson
    ! ---

    !n_states = nOO + 5
    !n_states_diag = n_states + 4
    !allocate(Om(nOO+nVV), R(nOO+nVV,n_states_diag))

    !supp_data_int = 1
    !allocate(supp_data_int(supp_data_int_size))
    !supp_data_int(1) = nS

    !supp_data_dbl_size = nS + nOrb*nOrb*nS + 1
    !allocate(supp_data_dbl(supp_data_dbl_size))
    !! scalars
    !supp_data_dbl(1) = eta
    !i_data = 1
    !! rho_RPA
    !do q = 1, nOrb
    !  do p = 1, nOrb
    !    do m = 1, nS
    !      i_data = i_data + 1
    !      supp_data_dbl(i_data) = rho_RPA(p,q,m)
    !    enddo
    !  enddo
    !enddo
    !! OmRPA
    !do m = 1, nS
    !  i_data = i_data + 1
    !  supp_data_dbl(i_data) = OmRPA(m)
    !enddo

    !call ppLR_davidson(ispin, TDA, nC, nO, nR, nOrb, nOO, nVV, &
    !                   1.d0,                                   & ! lambda
    !                   eGW(1),                                 &
    !                   0.d0,                                   & ! eF
    !                   ERI(1,1,1,1),                           &
    !                   supp_data_int(1), supp_data_int_size,   &
    !                   supp_data_dbl(1), supp_data_dbl_size,   &
    !                   Om(1), R(1,1), n_states, n_states_diag, "GW")

    !deallocate(Om, R)
    !deallocate(supp_data_dbl, supp_data_int)
    !stop

    ! ---

    call ppLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    !----------------------------------------------------!
    ! Compute the dynamical screening at the ppBSE level !
    !----------------------------------------------------!

    if(dBSE) &
        call RGW_ppBSE_dynamic_perturbation(ispin,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           Om1,X1,Y1,Om2,X2,Y2,KB_sta,KC_sta,KD_sta)


    deallocate(Om1,X1,Y1,Om2,X2,Y2)
  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    write(*,*) '****************'
    write(*,*) '*** Triplets ***'
    write(*,*) '****************'
    write(*,*) 

    ispin = 2
    EcBSE(ispin) = 0d0

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV))
    allocate(Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    ! Compute BSE excitation energies

    ! ---
    ! LAPACK
    ! ---

    allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
    allocate(KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

    call RGW_ppBSE_static_kernel_C(ispin,eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,OmRPA,rho_RPA,KC_sta)
    call RGW_ppBSE_static_kernel_D(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,OmRPA,rho_RPA,KD_sta)
    if(.not.TDA) call RGW_ppBSE_static_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,KB_sta)

    call ppLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
    call ppLR_D(ispin,nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
    if(.not.TDA) call ppLR_B(ispin,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE(ispin))
    deallocate(Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)

    !print*, 'LAPACK:'
    !print*, Om2
    !print*, Om1

    ! ---
    ! Davidson
    ! ---

    !n_states = nOO + 5
    !n_states_diag = n_states + 4
    !allocate(Om(nOO+nVV), R(nOO+nVV,n_states_diag))

    !supp_data_int = 1
    !allocate(supp_data_int(supp_data_int_size))
    !supp_data_int(1) = nS

    !supp_data_dbl_size = nS + nOrb*nOrb*nS + 1
    !allocate(supp_data_dbl(supp_data_dbl_size))
    !! scalars
    !supp_data_dbl(1) = eta
    !i_data = 1
    !! rho_RPA
    !do q = 1, nOrb
    !  do p = 1, nOrb
    !    do m = 1, nS
    !      i_data = i_data + 1
    !      supp_data_dbl(i_data) = rho_RPA(p,q,m)
    !    enddo
    !  enddo
    !enddo
    !! OmRPA
    !do m = 1, nS
    !  i_data = i_data + 1
    !  supp_data_dbl(i_data) = OmRPA(m)
    !enddo

    !call ppLR_davidson(ispin, TDA, nC, nO, nR, nOrb, nOO, nVV, &
    !                   1.d0,                                   & ! lambda
    !                   eGW(1),                                 &
    !                   0.d0,                                   & ! eF
    !                   ERI(1,1,1,1),                           &
    !                   supp_data_int(1), supp_data_int_size,   &
    !                   supp_data_dbl(1), supp_data_dbl_size,   &
    !                   Om(1), R(1,1), n_states, n_states_diag, "GW")

    !deallocate(Om, R)
    !deallocate(supp_data_dbl, supp_data_int)
    !stop

    ! ---

    EcBSE(ispin) = 3d0*EcBSE(ispin)

    call ppLR_transition_vectors(.false.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    !----------------------------------------------------!
    ! Compute the dynamical screening at the ppBSE level !
    !----------------------------------------------------!

    if(dBSE) &
        call RGW_ppBSE_dynamic_perturbation(ispin,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           Om1,X1,Y1,Om2,X2,Y2,KB_sta,KC_sta,KD_sta)

    deallocate(Om1,X1,Y1,Om2,X2,Y2)

  end if

end subroutine 
