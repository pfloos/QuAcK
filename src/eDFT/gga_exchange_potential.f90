subroutine gga_exchange_potential(DFA,nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,dAO,&
                                  rho,drho,Cx_choice,Fx)

! Select GGA exchange functional for potential calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: DFA
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nCC
  double precision,intent(in)   :: aCC(nCC,nEns-1)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)
  integer,intent(in)            :: Cx_choice

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Select GGA exchange functional

  select case (DFA)

    case (1)

      call G96_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)
      
    case (2)
    
      call B88_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

    case (3)
    
      call PBE_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

    case (4)

      call CC_B88_gga_exchange_potential(nEns,wEns,nCC,aCC,nGrid,weight,nBas,AO,dAO,rho,drho,&
                                         Cx_choice,Fx)

    case default

      call print_warning('!!! GGA exchange potential not available !!!')
      stop

  end select

end subroutine gga_exchange_potential
