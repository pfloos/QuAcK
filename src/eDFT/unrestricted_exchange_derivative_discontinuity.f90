subroutine unrestricted_exchange_derivative_discontinuity(rung,DFA,nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,drhow,ExDD)

! Compute the exchange part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  character(len=12),intent(in)  :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)

! Local variables


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

  select case (rung)

!   Hartree calculation

    case(0) 

      ExDD(:) = 0d0

!   LDA functionals

    case(1) 

      call unrestricted_lda_exchange_derivative_discontinuity(DFA,nEns,wEns(:),aCC_w1,aCC_w2,nGrid,weight(:),rhow(:),ExDD(:))

!   GGA functionals

    case(2) 

      call unrestricted_gga_exchange_derivative_discontinuity(DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),ExDD(:))

!   Hybrid functionals

    case(4) 

      call print_warning('!!! exchange part of derivative discontinuity NYI for hybrids !!!')
      stop

!   Hartree-Fock calculation

    case(666) 

      ExDD(:) = 0d0

  end select
 
end subroutine unrestricted_exchange_derivative_discontinuity
