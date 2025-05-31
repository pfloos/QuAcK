subroutine complex_qsRGW_plot_self_energy(nBas,eta,nC,nO,nV,nR,nS,eGW)

! Dump several GW quantities for external plotting

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: eGW(nBas)

! Local variables

  integer                       :: p,g
  integer                       :: nP
  character(len=256)            :: fmtP
  integer                       :: nGrid
  double precision              :: wmin,wmax,dw
  double precision,allocatable  :: w(:)
  double precision,allocatable  :: A(:,:)

! Construct grid

  nGrid = 20000
  allocate(w(nGrid),A(nBas,nGrid))

! Minimum and maximum frequency values

  wmin = -2d0
  wmax = +2d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  end do


! Compute spectral function
  do g=1,nGrid
    do p=nC+1,nO
      A(p,g) = abs(0d0 - aimag(eGW(p))-eta)&
              /((w(g) - real(eGW(p)))**2 + (0d0 - aimag(eGW(p))- eta)**2)
    end do
  end do
  do g=1,nGrid
    do p=nO+1,nBas-nR
      A(p,g) = abs(0d0 - aimag(eGW(p))+eta)&
              /((w(g) - real(eGW(p)))**2 + (0d0 - aimag(eGW(p))+eta)**2)
    end do
  end do

  A(:,:) = A(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=11 ,file='qsRGW_A.dat')

  nP = nBas - nR - nC 
  write(fmtP, '(A,I0,A)') '(F12.6,1X,', nP, '(F12.6,1X))'
  do g=1,nGrid
    write(11,fmtP) w(g)*HaToeV,(A(p,g),p=nC+1,nBas-nR)
  end do

! Closing files and deallocation
  close(unit=11)
  deallocate(w,A)

end subroutine 
