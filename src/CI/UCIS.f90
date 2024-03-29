subroutine UCIS(dotest,spin_conserved,spin_flip,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb, &
                dipole_int_aa,dipole_int_bb,eHF,cHF,S)

! Perform configuration interaction single calculation`

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart,nspin)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart,nspin)

! Local variables

  logical                       :: dump_matrix = .false.
  logical                       :: dump_trans  = .false.
  integer                       :: ispin
  double precision              :: lambda

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: A_sc(:,:)
  double precision,allocatable  :: Om_sc(:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: A_sf(:,:)
  double precision,allocatable  :: Om_sf(:)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Unrestrictd CIS Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Adiabatic connection scaling

  lambda = 1d0

!----------------------------!
! Spin-conserved transitions !
!----------------------------!

  if(spin_conserved) then

    ispin = 1

    ! Memory allocation

    nS_aa = nS(1)
    nS_bb = nS(2)
    nS_sc = nS_aa + nS_bb

    allocate(A_sc(nS_sc,nS_sc),Om_sc(nS_sc))

    call phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A_sc)
 
    if(dump_matrix) then
      print*,'CIS matrix (spin-conserved transitions)'
      call matout(nS_sc,nS_sc,A_sc)
      write(*,*)
    end if

    call diagonalize_matrix(nS_sc,A_sc,Om_sc)
    A_sc(:,:) = transpose(A_sc)
    call print_excitation_energies('CIS@UHF','spin-conserved',nS_sc,Om_sc)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc,dipole_int_aa,dipole_int_bb, &
                                  cHF,S,Om_sc,A_sc,A_sc)
 
    if(dump_trans) then
      print*,'Spin-conserved CIS transition vectors'
      call matout(nS_sc,nS_sc,A_sc)
      write(*,*)
    end if

    ! Testing zone

    if(dotest) then

      call dump_test_value('U','CIS singlet excitation energy',Om_sc(1))

    end if

    deallocate(A_sc,Om_sc)

  end if

!-----------------------!
! Spin-flip transitions !
!-----------------------!

  if(spin_flip) then

    ispin = 2

    ! Memory allocation

    nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
    nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
    nS_sf = nS_ab + nS_ba

    allocate(A_sf(nS_sf,nS_sf),Om_sf(nS_sf))
    
    call phULR_A(ispin,.false.,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,lambda,eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,A_sf)

    if(dump_matrix) then
      print*,'CIS matrix (spin-conserved transitions)'
      call matout(nS_sf,nS_sf,A_sf)
      write(*,*)
    end if

    call diagonalize_matrix(nS_sf,A_sf,Om_sf)
    A_sf(:,:) = transpose(A_sf)
    call print_excitation_energies('CIS@UHF','spin-flip',nS_sf,Om_sf)
    call phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf,dipole_int_aa,dipole_int_bb, &
                                  cHF,S,Om_sf,A_sf,A_sf)
 
    if(dump_trans) then
      print*,'Spin-flip CIS transition vectors'
      call matout(nS_sf,nS_sf,A_sf)
      write(*,*)
    end if

    ! Testing zone

    if(dotest) then

      call dump_test_value('U','CIS triplet excitation energy',Om_sf(1))

    end if

    deallocate(A_sf,Om_sf)

  end if

end subroutine 
