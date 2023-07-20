subroutine GTpp_phBSE(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa,   &
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

  logical                       :: dRPA

  integer                       :: ispin
  integer                       :: iblock

  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: TAab(:,:),TBab(:,:)
  double precision,allocatable  :: TAaa(:,:),TBaa(:,:)
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),TAab(nS,nS),TBab(nS,nS),TAaa(nS,nS),TBaa(nS,nS), & 
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!---------------------------------------!
! Compute T-matrix for alpha-beta block !
!---------------------------------------!

  ispin  = 1
  iblock = 3

  allocate(Bpp(nVVab,nOOab),Cpp(nVVab,nVVab),Dpp(nOOab,nOOab))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVab,1d0,eT,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOab,1d0,eT,ERI,Dpp)

  call ppLR(TDA_T,nOOab,nVVab,Bpp,Cpp,Dpp,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

               call GTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TAab)
  if(.not.TDA) call GTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,1d0,Om1ab,rho1ab,Om2ab,rho2ab,TBab)

!----------------------------------------!
! Compute T-matrix for alpha-alpha block !
!----------------------------------------!

  ispin  = 2
  iblock = 4

  allocate(Bpp(nVVaa,nOOaa),Cpp(nVVaa,nVVaa),Dpp(nOOaa,nOOaa))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVaa,1d0,eT,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOaa,1d0,eT,ERI,Dpp)

  call ppLR(TDA_T,nOOaa,nVVaa,Bpp,Cpp,Dpp,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

               call GTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOaa,nVVaa,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TAaa)
  if(.not.TDA) call GTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOaa,nVVaa,1d0,Om1aa,rho1aa,Om2aa,rho2aa,TBaa)

!------------------!
! Singlet manifold !
!------------------!

 if(singlet) then

    ispin = 1

    ! Compute BSE singlet excitation energies

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + TAaa(:,:) + TAab(:,:) 
    if(.not.TDA) Bph(:,:) = Bph(:,:) + TBaa(:,:) + TBab(:,:) 

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('phBSE@GTpp  ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    if(dBSE) then
 
      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        print*, ' Iterative dynamical correction for BSE@GT NYI'

      else
 
        call GTpp_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa, &
                                             Om1ab,Om2ab,Om1aa,Om2aa,rho1ab,rho2ab,rho1aa,rho2aa,eT,eGT, & 
                                             dipole_int,OmBSE,XpY_BSE,XmY_BSE,TAab,TAaa)
      end if
 
    end if

  end if

!------------------!
! Triplet manifold !
!------------------!

 if(triplet) then

    ispin = 2

    ! Compute BSE triplet excitation energies

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + TAaa(:,:) - TAab(:,:) 
    if(.not.TDA) Bph(:,:) = Bph(:,:) + TBaa(:,:) - TBab(:,:) 

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('phBSE@GTpp  ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    if(dBSE) then
 
      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        print*, ' Iterative dynamical correction for BSE@GT NYI'

      else
 
        call GTpp_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa, &
                                             Om1ab,Om2ab,Om1aa,Om2aa,rho1ab,rho2ab,rho1aa,rho2aa,eT,eGT, & 
                                             dipole_int,OmBSE,XpY_BSE,XmY_BSE,TAab,TAaa)
      end if
 
    end if

  end if

end subroutine 
