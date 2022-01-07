subroutine Bethe_Salpeter_Tmatrix(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                                  Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,rho1s,rho2s,Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,rho1t,rho2t, & 
                                  ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

  double precision,intent(in)   :: Omega1s(nVVs)
  double precision,intent(in)   :: X1s(nVVs,nVVs)
  double precision,intent(in)   :: Y1s(nOOs,nVVs)
  double precision,intent(in)   :: Omega2s(nOOs)
  double precision,intent(in)   :: X2s(nVVs,nOOs)
  double precision,intent(in)   :: Y2s(nOOs,nOOs)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs)
  double precision,intent(in)   :: Omega1t(nVVt)
  double precision,intent(in)   :: X1t(nVVt,nVVt)
  double precision,intent(in)   :: Y1t(nOOt,nVVt)
  double precision,intent(in)   :: Omega2t(nOOt)
  double precision,intent(in)   :: X2t(nVVt,nOOt)
  double precision,intent(in)   :: Y2t(nOOt,nOOt) 
  double precision,intent(in)   :: rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: rho2t(nBas,nBas,nOOt)

! Local variables

  integer                       :: ispin
  integer                       :: iblock

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: TA(:,:),TB(:,:)
  double precision,allocatable  :: OmBSE(:,:)
  double precision,allocatable  :: XpY_BSE(:,:,:)
  double precision,allocatable  :: XmY_BSE(:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(TA(nS,nS),TB(nS,nS),OmBSE(nS,nspin),XpY_BSE(nS,nS,nspin),XmY_BSE(nS,nS,nspin))

! Initialize T matrix

  TA(:,:) = 0d0
  TB(:,:) = 0d0

!----------------------------------------------
! Compute T-matrix for alpha-beta block
!----------------------------------------------

  ispin  = 1
  iblock = 3

  call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eT,ERI,  &
                          Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,EcRPA(ispin))

! call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

               call static_Tmatrix_A(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,ERI,Omega1s,rho1s,Omega2s,rho2s,TA)
  if(.not.TDA) call static_Tmatrix_B(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,ERI,Omega1s,rho1s,Omega2s,rho2s,TB)

  print*,'aa block of TA'
  call matout(nS,nS,TA)
  print*,'aa block of TB'
  call matout(nS,nS,TB)

!----------------------------------------------
! Compute T-matrix for alpha-alpha block
!----------------------------------------------

  ispin  = 2
  iblock = 4

  call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eT,ERI,  &
                          Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,EcRPA(ispin))

! call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

               call static_Tmatrix_A(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,ERI,Omega1t,rho1t,Omega2t,rho2t,TA)
  if(.not.TDA) call static_Tmatrix_B(ispin,eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,ERI,Omega1t,rho1t,Omega2t,rho2t,TB)

  print*,'aa+ab block of TA'
  call matout(nS,nS,TA)
  print*,'aa+ab block of TB'
  call matout(nS,nS,TB)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE singlet excitation energies

    call linear_response_Tmatrix(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TA,TB, &
                                 EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then

        print*, ' Iterative dynamical correction for BSE@GT NYI'
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else

        call Bethe_Salpeter_Tmatrix_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,Omega1s,Omega2s,rho1s,rho2s, &
                                                         eT,eGT,dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      end if

    end if

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE triplet excitation energies

    call linear_response_Tmatrix(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TA,TB, &
                                 EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then

    
        print*, ' Iterative dynamical correction for BSE@GT NYI'
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else

        call Bethe_Salpeter_Tmatrix_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,Omega1t,Omega2t,rho1t,rho2t, &
                                                         eT,eGT,dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      end if

    end if

  end if

end subroutine Bethe_Salpeter_Tmatrix
