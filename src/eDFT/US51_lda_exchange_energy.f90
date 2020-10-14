subroutine US51_lda_exchange_energy(nGrid,weight,rho,ExLDA)

  use xc_f90_lib_m

! Compute Slater's LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer(8)                    :: nGri8 
  integer                       :: iG
  double precision              :: r
  double precision,allocatable  :: Ex(:)

  TYPE(xc_f90_func_t) :: xc_func
  TYPE(xc_f90_func_info_t) :: xc_info
  integer :: func_id = 1

! Output variables

  double precision              :: ExLDA

! Memory allocation
  
  nGri8 = int(nGrid,8)
  print*,nGri8
  allocate(Ex(nGrid))

  call xc_f90_func_init(xc_func, func_id, XC_POLARIZED)
  xc_info = xc_f90_func_get_info(xc_func)
  call xc_f90_lda_exc(xc_func, nGri8, rho(1), Ex(1))

  ExLDA = 0d0

! do iG=1,nGrid

! write(*,"(F8.6,1X,F9.6)") rho(iG), Ex(iG)

!     ExLDA = ExLDA  + weight(iG)*Ex(iG)

! enddo

  call xc_f90_func_end(xc_func)


end subroutine US51_lda_exchange_energy
