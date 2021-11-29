subroutine unrestricted_exchange_derivative_discontinuity(rung,DFA,nEns,wEns,nCC,aCC,nGrid,weight,rhow,drhow,&
                                                          Cx_choice,doNcentered,kappa,ExDD)

! Compute the exchange part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: rung
  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  integer,intent(in)            :: Cx_choice
  logical,intent(in)            :: doNcentered
  double precision,intent(in)   :: kappa(nEns)

! Local variables


! Output variables

  double precision,intent(out)  :: ExDD(nEns)

  select case (rung)

!   Hartree calculation

    case(0) 

      ExDD(:) = 0d0

!   LDA functionals

    case(1) 

      call unrestricted_lda_exchange_derivative_discontinuity(DFA,nEns,wEns(:),nCC,aCC,nGrid,weight(:),&
                                                              rhow(:),Cx_choice,doNcentered,kappa,ExDD(:))
!   GGA functionals

    case(2) 

      call unrestricted_gga_exchange_derivative_discontinuity(DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),ExDD(:))

!   MGGA functionals

    case(3) 

      call unrestricted_mgga_exchange_derivative_discontinuity(DFA,nEns,wEns(:),nGrid,weight(:),rhow(:),drhow(:,:),ExDD(:))

!   Hybrid functionals

    case(4) 

      call unrestricted_hybrid_exchange_derivative_discontinuity(DFA,nEns,wEns(:),nCC,aCC,nGrid,weight(:),&
                                                                 rhow(:),Cx_choice,doNcentered,kappa,ExDD(:))

  end select
 
end subroutine unrestricted_exchange_derivative_discontinuity
