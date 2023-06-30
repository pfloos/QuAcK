double precision function dSigmaC(x,w,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,regularize)

! Compute the derivative of the correlation part of the self-energy

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
  logical,intent(in)            :: regularize

! Local variables

  integer                       :: i,j,a,b,p,jb
  double precision              :: eps
  double precision              :: Dpijb,Dpajb

! Initialize 

  dSigmaC = 0d0

  if (regularize) then
  ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        do jb=1,nS
           eps = w - e(i) + Omega(jb)
           Dpijb = e(p) - e(i) + Omega(jb)
           dSigmaC = dSigmaC - 2d0*rho(p,i,jb)**2*(1d0-exp(-2*eta*Dpijb*Dpijb))/(eps**2)
        enddo
     enddo
  ! Virtual part of the correlation self-energy
     do a=nO+1,nBas-nR
        do jb=1,nS
           eps = w - e(a) - Omega(jb)
           Dpajb = e(p) - e(a) - Omega(jb)
           dSigmaC = dSigmaC - 2d0*rho(p,a,jb)**2*(1d0-exp(-2*eta*Dpajb*Dpajb))/(eps**2)
        enddo
     enddo

  else
   ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
           do b=nO+1,nBas-nR
              jb = jb + 1
              eps = w - e(i) + Omega(jb)
              dSigmaC = dSigmaC - 2d0*rho(x,i,jb)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
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
              dSigmaC = dSigmaC - 2d0*rho(x,a,jb)**2*(eps**2 - eta**2)/(eps**2 + eta**2)**2
           enddo
        enddo
     enddo
  end if

end function dSigmaC
