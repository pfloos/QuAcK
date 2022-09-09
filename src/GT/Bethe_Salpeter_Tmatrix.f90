subroutine Bethe_Salpeter_Tmatrix(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa,   &
                                  Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, & 
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

  integer,intent(in)            :: nOOab
  integer,intent(in)            :: nOOaa
  integer,intent(in)            :: nVVab
  integer,intent(in)            :: nVVaa

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

  double precision,intent(in)   :: Om1ab(nVVab)
  double precision,intent(in)   :: X1ab(nVVab,nVVab)
  double precision,intent(in)   :: Y1ab(nOOab,nVVab)
  double precision,intent(in)   :: Om2ab(nOOab)
  double precision,intent(in)   :: X2ab(nVVab,nOOab)
  double precision,intent(in)   :: Y2ab(nOOab,nOOab)
  double precision,intent(in)   :: rho1ab(nBas,nBas,nVVab)
  double precision,intent(in)   :: rho2ab(nBas,nBas,nOOab)
  double precision,intent(in)   :: Om1aa(nVVaa)
  double precision,intent(in)   :: X1aa(nVVaa,nVVaa)
  double precision,intent(in)   :: Y1aa(nOOaa,nVVaa)
  double precision,intent(in)   :: Om2aa(nOOaa)
  double precision,intent(in)   :: X2aa(nVVaa,nOOaa)
  double precision,intent(in)   :: Y2aa(nOOaa,nOOaa) 
  double precision,intent(in)   :: rho1aa(nBas,nBas,nVVaa)
  double precision,intent(in)   :: rho2aa(nBas,nBas,nOOaa)

! Local variables

  integer                       :: ispin
  integer                       :: iblock

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: TAab(:,:),TBab(:,:)
  double precision,allocatable  :: TAaa(:,:),TBaa(:,:)
  double precision,allocatable  :: OmBSE(:,:)
  double precision,allocatable  :: XpY_BSE(:,:,:)
  double precision,allocatable  :: XmY_BSE(:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(TAab(nS,nS),TBab(nS,nS),TAaa(nS,nS),TBaa(nS,nS), & 
           OmBSE(nS,nspin),XpY_BSE(nS,nS,nspin),XmY_BSE(nS,nS,nspin))

!---------------------------------------!
! Compute T-matrix for alpha-beta block !
!---------------------------------------!

  ispin  = 1
  iblock = 3

  call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,eT,ERI,  &
                          Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))

               call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TAab)
  if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TBab)

!----------------------------------------!
! Compute T-matrix for alpha-alpha block !
!----------------------------------------!

  ispin  = 2
  iblock = 4

  call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,eT,ERI,  &
                          Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

               call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOaa,nVVaa,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TAaa)
  if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOaa,nVVaa,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TBaa)

!------------------!
! Singlet manifold !
!------------------!

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE singlet excitation energies

    call linear_response_BSE(ispin,.false.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TAab+TAaa,TBab+TBaa, &
                             EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    if(dBSE) then
 
      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        print*, ' Iterative dynamical correction for BSE@GT NYI'
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else
 
        call Bethe_Salpeter_Tmatrix_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa, &
                                                         Om1ab,Om2ab,Om1aa,Om2aa,rho1ab,rho2ab,rho1aa,rho2aa,eT,eGT, & 
                                                         dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),TAab,TAaa)
      end if
 
    end if

  end if

!------------------!
! Triplet manifold !
!------------------!

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE triplet excitation energies

    call linear_response_BSE(ispin,.false.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TAaa-TAab,TBaa-TBab, &
                             EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    if(dBSE) then
 
      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        print*, ' Iterative dynamical correction for BSE@GT NYI'
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else
 
        call Bethe_Salpeter_Tmatrix_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa, &
                                                         Om1ab,Om2ab,Om1aa,Om2aa,rho1ab,rho2ab,rho1aa,rho2aa,eT,eGT, & 
                                                         dipole_int,OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),TAab,TAaa)
      end if
 
    end if

  end if

end subroutine Bethe_Salpeter_Tmatrix
