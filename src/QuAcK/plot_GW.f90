subroutine plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eGW,Omega,rho,rhox)

! Dump several GW quantities for external plotting

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas),eGW(nBas),Omega(nS),rho(nBas,nBas,nS),rhox(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,x,jb,g
  integer                       :: nGrid
  double precision              :: eps,eta,wmin,wmax,dw
  double precision,allocatable  :: w(:),SigC(:,:),Z(:,:),S(:,:)

! Infinitesimal

  eta = 1d-3

! Construct grid

  nGrid = 1000
  allocate(w(nGrid),SigC(nBas,nGrid),Z(nBas,nGrid),S(nBas,nGrid))

! Initialize 

  SigC(:,:) = 0d0
  Z(:,:)    = 0d0

! Minimum and maximum frequency values

  wmin = -5d0
  wmax = +5d0
  dw = (wmax - wmin)/dble(ngrid)

  do g=1,nGrid
    w(g) = wmin + dble(g)*dw
  enddo

! Occupied part of the self-energy and renormalization factor
 
  do g=1,nGrid
    do x=nC+1,nBas-nR
      do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = w(g) - eHF(i) + Omega(jb)
            SigC(x,g) = SigC(x,g) + 2d0*rho(x,i,jb)**2*eps/(eps**2 + eta**2)
            Z(x,g) = Z(x,g) + 2d0*rho(x,i,jb)**2/eps**2
          enddo
        enddo
      enddo
    enddo
  enddo
 
! Virtual part of the self-energy and renormalization factor
 
  do g=1,nGrid
    do x=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = w(g) - eHF(a) - Omega(jb)
            SigC(x,g) = SigC(x,g) + 2d0*rho(x,a,jb)**2*eps/(eps**2 + eta**2)
            Z(x,g) = Z(x,g) + 2d0*rho(x,a,jb)**2/eps**2
          enddo
        enddo
      enddo
    enddo
  enddo

  Z(:,:) = 1d0/(1d0 + Z(:,:))

! Compute spectral function

  do g=1,nGrid
    do x=nC+1,nBas-nR
      S(x,g) = eta/((w(g) - eHF(x) - SigC(x,g))**2 + eta**2)
    enddo
  enddo

  S(:,:) = S(:,:)/pi

! Dump quantities in files as a function of w

  open(unit=8 ,file='plot/grid.dat')
  open(unit=9 ,file='plot/SigC.dat')
  open(unit=10 ,file='plot/Z.dat')
  open(unit=11 ,file='plot/A.dat')

  do g=1,nGrid
    write(8 ,*) w(g)*HaToeV,(SigC(x,g)*HaToeV,x=1,nBas)
    write(9 ,*) w(g)*HaToeV,((w(g)-eHF(x))*HaToeV,x=1,nBas)
    write(10,*) w(g)*HaToeV,(Z(x,g),x=1,nBas)
    write(11,*) w(g)*HaToeV,(S(x,g),x=1,nBas)
  enddo

! Closing files

  close(unit=8)
  close(unit=9)
  close(unit=10)
  close(unit=11)

end subroutine plot_GW
