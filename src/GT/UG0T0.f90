subroutine UG0T0(doACFDT,exchange_kernel,doXBS,BSE,TDA_T,TDA,dBSE,dTDA,evDyn, &
                 spin_conserved,spin_flip,linearize,eta,regularize,nBas,nC,nO,nV, &
                 nR,nS,ENuc,EUHF,ERI,ERI_aaaa,ERI_aabb,ERI_bbbb, &
                 dipole_int_aa,dipole_int_bb,PHF,cHF,eHF,Vxc,eG0T0)

! Perform one-shot calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin) 
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: Vxc(nBas,nspin)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas) 
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin,is
  integer                       :: iblock 
  integer                       :: nH_sc,nH_sf,nHaa,nHab,nHbb
  integer                       :: nP_sc,nP_sf,nPaa,nPab,nPbb
  double precision              :: EcRPA(nspin),Ecaa,Ecbb
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Omega1ab(:),Omega1aa(:),Omega1bb(:)
  double precision,allocatable  :: X1ab(:,:),X1aa(:,:),X1bb(:,:)
  double precision,allocatable  :: Y1ab(:,:),Y1aa(:,:),Y1bb(:,:)
  double precision,allocatable  :: rho1ab(:,:,:),rho1aa(:,:,:),rho1bb(:,:,:)
  double precision,allocatable  :: Omega2ab(:),Omega2aa(:),Omega2bb(:)
  double precision,allocatable  :: X2ab(:,:),X2aa(:,:),X2bb(:,:)
  double precision,allocatable  :: Y2ab(:,:),Y2aa(:,:),Y2bb(:,:)
  double precision,allocatable  :: rho2ab(:,:,:),rho2aa(:,:,:),rho2bb(:,:,:)
  double precision,allocatable  :: SigX(:,:)
  double precision,allocatable  :: SigT(:,:)
  double precision,allocatable  :: Z(:,:)

! Output variables

  double precision,intent(out)  :: eG0T0(nBas,nspin)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0 calculation           |'
  write(*,*)'|         *** Unrestricted version ***         |'
  write(*,*)'************************************************'
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

! Memory allocation

  allocate(Omega1ab(nPab),X1ab(nPab,nPab),Y1ab(nHab,nPab), & 
           Omega2ab(nHab),X2ab(nPab,nHab),Y2ab(nHab,nHab), & 
           rho1ab(nBas,nBas,nPab),rho2ab(nBas,nBas,nHab), & 
           Omega1aa(nPaa),X1aa(nPaa,nPaa),Y1aa(nHaa,nPaa), & 
           Omega2aa(nHaa),X2aa(nPaa,nHaa),Y2aa(nHaa,nHaa), & 
           rho1aa(nBas,nBas,nPaa),rho2aa(nBas,nBas,nHaa), & 
           Omega1bb(nPbb),X1bb(nPbb,nPbb),Y1bb(nHbb,nPbb), &
           Omega2bb(nPbb),X2bb(nPbb,nPbb),Y2bb(nHbb,nPbb), &
           rho1bb(nBas,nBas,nPbb),rho2bb(nBas,nBas,nHbb), &
           SigX(nBas,nspin),SigT(nBas,nspin),Z(nBas,nspin))

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  ispin  = 1
  iblock = 3
! iblock = 1

! Compute linear response

  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPab,nHaa,nHab,nHbb,nHab,1d0,eHF,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1ab,X1ab,Y1ab, &
                                       Omega2ab,X2ab,Y2ab,EcRPA(ispin)) 
call matout(nHab,nPab,Y1ab)  
! EcRPA(ispin) = 1d0*EcRPA(ispin)

  call print_excitation('pp-RPA (N+2)',iblock,nPab,Omega1ab(:))
  call print_excitation('pp-RPA (N-2)',iblock,nHab,Omega2ab(:))
 
!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

  ispin  = 2
  iblock = 4

! Compute linear response

  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPaa,nHaa,nHab,nHbb,nHaa,1d0,eHF,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1aa,X1aa,Y1aa, &
                                       Omega2aa,X2aa,Y2aa,EcRPA(ispin))  
  
! EcRPA(ispin) = 2d0*EcRPA(ispin)
! EcRPA(ispin) = 3d0*EcRPA(ispin)

  call print_excitation('pp-RPA (N+2)',iblock,nPaa,Omega1aa(:))
  call print_excitation('pp-RPA (N-2)',iblock,nHaa,Omega2aa(:))

!----------------------------------------------
! beta-beta block
!----------------------------------------------

  ispin  = 2
  iblock = 7

! Compute linear response

  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPbb,nHaa,nHab,nHbb,nHbb,1d0,eHF,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1bb,X1bb,Y1bb, &
                                       Omega2bb,X2bb,Y2bb,EcRPA(ispin))

! EcRPA(ispin) = 2d0*EcRPA(ispin)
! EcRPA(ispin) = 3d0*EcRPA(ispin)

  call print_excitation('pp-RPA (N+2)',iblock,nPbb,Omega1bb(:))
  call print_excitation('pp-RPA (N-2)',iblock,nHbb,Omega2bb(:))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  EcGM    = 0d0
  SigT(:,:) = 0d0
  Z(:,:)    = 0d0

!alpha-beta block
  
  iblock = 3

  call unrestricted_excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nHab,nPab, &
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,X1ab,Y1ab, &
                                               rho1ab,X2ab,Y2ab,rho2ab)
!alpha-alpha block

  iblock = 4
  
  call unrestricted_excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nHaa,nPaa, &
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,X1aa,Y1aa, &
                                               rho1aa,X2aa,Y2aa,rho2aa)

!beta-beta block 
  
  iblock = 7

  call unrestricted_excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nHbb,nPbb, &
                                               ERI_aaaa,ERI_aabb,ERI_bbbb,X1bb,Y1bb, &
                                               rho1bb,X2bb,Y2bb,rho2bb)

  call unrestricted_self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,nPaa,&
                                             nPab,nPbb,eHF,Omega1aa,Omega1ab,Omega1bb,&
                                             rho1aa,rho1ab,rho1bb,Omega2aa,Omega2ab,&
                                             Omega2bb,rho2aa,rho2ab,rho2bb,EcGM,SigT)

  call unrestricted_renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nHaa,nHab,nHbb,&
                                                   nPaa,nPab,nPbb,eHF,Omega1aa,Omega1ab,&
                                                   Omega1bb,rho1aa,rho1ab,rho1bb, &
                                                   Omega2aa,Omega2ab,Omega2bb,rho2aa, &
                                                   rho2ab,rho2bb,Z) 


  Z(:,:) = 1d0/(1d0 - Z(:,:))

!----------------------------------------------
! Compute the exchange part of the self-energy
!----------------------------------------------
  
  do is=1,nspin
    call self_energy_exchange_diag(nBas,cHF(:,:,is),PHF(:,:,is),ERI,SigX(:,is))
  end do 
!call matout(nBas,nspin,SigX)
!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------

  if(linearize) then

!   eG0T0(:) = eHF(:) + Z(:)*SigT(:)
    eG0T0(:,:) = eHF(:,:) + Z(:,:)*(SigX(:,:) + SigT(:,:) - Vxc(:,:))
    
!    call matout(nBas,1,SigX)
!    call matout(nBas,1,Vxc)
!    call matout(nBas,1,eG0T0(:,1)*HaToeV)
!    call matout(nBas,nspin,SigT*HaToeV) 
  else
  
    eG0T0(:,:) = eHF(:,:) + SigX(:,:) + SigT(:,:) - Vxc(:,:)

  end if

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

!alpha-beta block

  ispin  = 1
  iblock = 3 

  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPab,nHaa,nHab,nHbb,nHab,1d0,eG0T0,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1ab,X1ab,Y1ab, &
                                       Omega2ab,X2ab,Y2ab,EcRPA(ispin))

!alpha-alpha block
 
  ispin  = 2
  iblock = 4 
  
  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPaa,nHaa,nHab,nHbb,nHaa,1d0,eG0T0,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1aa,X1aa,Y1aa, &
                                       Omega2aa,X2aa,Y2aa,EcRPA(ispin)) 
 
  Ecaa = EcRPA(2)

!beta-beta block
 
  iblock = 7
  
  call unrestricted_linear_response_pp(iblock,TDA,nBas,nC,nO,nV,nR,nPaa,nPab,nPbb, &
                                       nPbb,nHaa,nHab,nHbb,nHbb,1d0,eG0T0,ERI_aaaa, &
                                       ERI_aabb,ERI_bbbb,Omega1bb,X1bb,Y1bb, &
                                       Omega2bb,X2bb,Y2bb,EcRPA(ispin))
  
  Ecbb = EcRPA(2) 
  EcRPA(2) = Ecaa + Ecbb
  EcRPA(1) = EcRPA(1) - EcRPA(2)
  EcRPA(2) = 3d0*EcRPA(2)

  call print_UG0T0(nBas,nO,eHF,ENuc,EUHF,SigT,Z,eG0T0,EcGM,EcRPA)

! Free memory

  deallocate(Omega1ab,X1ab,Y1ab,Omega2ab,X2ab,Y2ab,rho1ab,rho2ab, &
             Omega1aa,X1aa,Y1aa,Omega2aa,X2aa,Y2aa,rho1aa,rho2aa, &
             Omega1bb,X1bb,Y1bb,Omega2bb,X2bb,Y2bb,rho1bb,rho2bb)

end subroutine UG0T0
