double precision function GW_SigC(p,w,eta,nBas,nC,nO,nV,nR,nS,e,Om,rho,regularize)

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
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  logical,intent(in)            :: regularize

! Local variables

  integer                       :: i,a,jb
  double precision              :: eps
  double precision              :: Dpijb,Dpajb

! Initialize 

  GW_SigC = 0d0

  if (regularize) then
  ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        do jb=1,nS
           eps = w - e(i) + Om(jb)
           Dpijb = e(p) - e(i) + Om(jb)
           GW_SigC = GW_SigC + 2d0*rho(p,i,jb)**2*(1d0-exp(-2d0*eta*Dpijb*Dpijb))/eps
        enddo
     enddo
  ! Virtual part of the correlation self-energy
     do a=nO+1,nBas-nR
        do jb=1,nS
           eps = w - e(a) - Om(jb)
           Dpajb = e(p) - e(a) - Om(jb)
           GW_SigC = GW_SigC + 2d0*rho(p,a,jb)**2*(1d0-exp(-2d0*eta*Dpajb*Dpajb))/eps
        enddo
     enddo

     ! We add the static SRG term in the self-energy directly
     ! do i=nC+1,nO
     !    do jb=1,nS
     !       Dpijb = e(p) - e(i) + Om(jb)
     !       SigmaC = SigmaC + 2d0*rho(p,i,jb)**2*(exp(-2d0*eta*Dpijb*Dpijb)/Dpijb)
     !    enddo
     ! enddo
     ! do a=nO+1,nBas-nR
     !    do jb=1,nS
     !       Dpajb = e(p) - e(a) - Om(jb)
     !       SigmaC = SigmaC + 2d0*rho(p,a,jb)**2*(exp(-2d0*eta*Dpajb*Dpajb)/Dpajb)
     !    enddo
     ! enddo
  
  else
   ! Occupied part of the correlation self-energy
     do i=nC+1,nO
        do jb=1,nS
           eps = w - e(i) + Om(jb)
           GW_SigC = GW_SigC + 2d0*rho(p,i,jb)**2*eps/(eps**2 + eta**2)
        enddo
     enddo
  ! Virtual part of the correlation self-energy
     do a=nO+1,nBas-nR
        do jb=1,nS
           eps = w - e(a) - Om(jb)
           GW_SigC = GW_SigC + 2d0*rho(p,a,jb)**2*eps/(eps**2 + eta**2)
        enddo
     enddo
  end if 

end function 
