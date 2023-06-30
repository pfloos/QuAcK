double precision function SigmaC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,regularize)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
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

  integer                       :: i,a,jb
  double precision              :: eps
  double precision              :: Dpijb,Dpajb

! Initialize 

  SigmaC = 0d0


  if (regularize) then
  ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        do jb=1,nS
           eps = w - e(i) + Omega(jb)
           Dpijb = e(p) - e(i) + Omega(jb)
           SigmaC = SigmaC + 2d0*rho(p,i,jb)**2*(1d0-exp(-2d0*eta*Dpijb*Dpijb))/eps
        enddo
     enddo
  ! Virtual part of the correlation self-energy
     do a=nO+1,nBas-nR
        do jb=1,nS
           eps = w - e(a) - Omega(jb)
           Dpajb = e(p) - e(a) - Omega(jb)
           SigmaC = SigmaC + 2d0*rho(p,a,jb)**2*(1d0-exp(-2d0*eta*Dpajb*Dpajb))/eps
        enddo
     enddo

     ! We add the static SRG term in the self-energy directly
     ! do i=nC+1,nO
     !    do jb=1,nS
     !       Dpijb = e(p) - e(i) + Omega(jb)
     !       SigmaC = SigmaC + 2d0*rho(p,i,jb)**2*(exp(-2d0*eta*Dpijb*Dpijb)/Dpijb)
     !    enddo
     ! enddo
     ! do a=nO+1,nBas-nR
     !    do jb=1,nS
     !       Dpajb = e(p) - e(a) - Omega(jb)
     !       SigmaC = SigmaC + 2d0*rho(p,a,jb)**2*(exp(-2d0*eta*Dpajb*Dpajb)/Dpajb)
     !    enddo
     ! enddo
  
  else
   ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        do jb=1,nS
           eps = w - e(i) + Omega(jb)
           SigmaC = SigmaC + 2d0*rho(p,i,jb)**2*eps/(eps**2 + eta**2)
        enddo
     enddo
  ! Virtual part of the correlation self-energy
     do a=nO+1,nBas-nR
        do jb=1,nS
           eps = w - e(a) - Omega(jb)
           SigmaC = SigmaC + 2d0*rho(p,a,jb)**2*eps/(eps**2 + eta**2)
        enddo
     enddo
  end if 

end function SigmaC
