subroutine UG0F2(dotest,BSE,TDA,dBSE,dTDA,spin_conserved,spin_flip,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EUHF, &
                 ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,eHF)

! Perform unrestricted G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
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
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas,nspin)

! Local variables

  integer                       :: is
  integer                       :: ispin
  double precision              :: Ec(nsp)
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc

  double precision,allocatable  :: eGFlin(:,:)
  double precision,allocatable  :: eGF(:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Unrestricted G0F2 Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),eGF(nBas,nspin),eGFlin(nBas,nspin))

!---------------------!
! Compute self-energy !
!---------------------!

  if(regularize) then

    call UGF2_reg_self_energy_diag(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eHF,SigC,Z)

  else

    call UGF2_self_energy_diag(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eHF,SigC,Z)

  end if

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGFlin(:,:) = eHF(:,:) + Z(:,:)*SigC(:,:)

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGF(:,:) = eGFlin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

    write(*,*) '!!! Graphical solution NYI for UG0F2 !!!'
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '

    eGF(:,:) = eGFlin(:,:)      
 
  end if

! Compute MP2 correlation energy

  call UMP2(.false.,nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EUHF,eGF,Ec)

! Dump results

  call print_UG0F2(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGF,Ec)

! Perform BSE calculation

  if(BSE) then

    print*,'!!! BSE2 NYI for UG0F2 !!!'

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('U','UG0F2 correlation energy',Ec)
    call dump_test_value('U','UG0F2 HOMOa energy',eGF(nO(1),1))
    call dump_test_value('U','UG0F2 LUMOa energy',eGF(nO(1)+1,1))
    call dump_test_value('U','UG0F2 HOMOa energy',eGF(nO(2),2))
    call dump_test_value('U','UG0F2 LUMOa energy',eGF(nO(2)+1,2))

  end if

end subroutine 
