subroutine hartree_individual_energy(rung,nBas,ERI,J,Pw,P,EJ)

! Compute the exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: J(nBas,nBas)
  double precision,intent(in)   :: Pw(nBas,nBas)
  double precision,intent(in)   :: P(nBas,nBas)

! Local variables

  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EJ

  select case (rung)

!   Hartree calculation

    case(0) 

      call hartree_coulomb(nBas,P(:,:),ERI(:,:,:,:),J(:,:))
      EJ = 0.5d0*trace_matrix(nBas,matmul(P(:,:),J(:,:)))

!   LDA functionals

    case(1) 

      call hartree_coulomb(nBas,Pw(:,:),ERI(:,:,:,:),J(:,:))
      EJ =       trace_matrix(nBas,matmul(P(:,:),J(:,:))) &
         - 0.5d0*trace_matrix(nBas,matmul(Pw(:,:),J(:,:)))

!   GGA functionals

    case(2) 

      call print_warning('!!! Hartee individual energies NYI for GGAs !!!')
      stop

!   Hybrid functionals

    case(4) 

      call print_warning('!!! Hartree individual energies NYI for Hybrids !!!')
      stop

!   Hartree-Fock calculation

    case(666) 

      call hartree_coulomb(nBas,P(:,:),ERI(:,:,:,:),J(:,:))
      EJ = 0.5d0*trace_matrix(nBas,matmul(P(:,:),J(:,:))) 

  end select
 
end subroutine hartree_individual_energy
