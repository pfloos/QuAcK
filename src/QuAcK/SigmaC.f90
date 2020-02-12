double precision function SigmaC(x,w,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: x
  double precision,intent(in)   :: w
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,p,jb
  double precision              :: eps

! Initialize 

  SigmaC = 0d0

! Occupied part of the correlation self-energy

  do i=nC+1,nO
    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        eps = w - e(i) + Omega(jb)
        SigmaC = SigmaC + 2d0*rho(x,i,jb)**2*eps/(eps**2 + eta**2)
      enddo
    enddo
  enddo

! Virtual part of the correlation self-energy

  do a=nO+1,nBas-nR
    jb = 0
    do j=nC+1,nO
      do b=nO+1,nBas-nR
        jb = jb + 1
        eps = w - e(a) - Omega(jb)
        SigmaC = SigmaC + 2d0*rho(x,a,jb)**2*eps/(eps**2 + eta**2)
      enddo
    enddo
  enddo

end function SigmaC
