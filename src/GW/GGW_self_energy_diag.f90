subroutine GGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,e,Om,rho,EcGM,Sig,Z,ERI)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor

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
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,a,b,p,m
  double precision              :: num,eps,dem1,dem2

! Output variables

  double precision,intent(out)  :: Sig(nBas)
  double precision,intent(out)  :: Z(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  Sig(:) = 0d0
  Z(:)   = 0d0

!----------------!
! GW self-energy !
!----------------!
 
! ! Occupied part of the correlation self-energy

!   do p=nC+1,nBas-nR
!     do i=nC+1,nO
!       do m=1,nS

!         eps = e(p) - e(i) + Om(m)
!         num = rho(p,i,m)**2
!         Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
!         Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

!       end do
!     end do
!   end do

! ! Virtual part of the correlation self-energy

!   do p=nC+1,nBas-nR
!     do a=nO+1,nBas-nR
!       do m=1,nS

!         eps = e(p) - e(a) - Om(m)
!         num = rho(p,a,m)**2
!         Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
!         Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

!       end do
!     end do
!   end do

!-----------------------------------------------!
! Testing another way to compute GT self-energy !
!-----------------------------------------------!

  do p=nC+1,nBas-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nBas-nR

              eps = e(p) + e(a) - e(i) - e(j)
              num = ERI(p,a,i,j)**2

              Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
              Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

           end do
           do a=nO+1,nBas-nR

              do m=1,nS
                 num = - ERI(p,i,j,a) * rho(i,a,m) * rho(j,p,m)
                 dem1 = e(p) + e(a) - e(i) - e(j)
                 dem2 = e(p) - e(j) + Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2 - num/dem1/dem2/dem2

                 num = - ERI(p,i,j,a) * rho(i,a,m) * rho(j,p,m)
                 dem1 = e(p) + e(a) - e(i) - e(j)
                 dem2 = e(a) - e(i) + Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
                 
                 num = - ERI(p,a,j,i) * rho(a,i,m) * rho(j,p,m)
                 dem1 = e(p) - e(j) + Om(m)
                 dem2 = e(a) - e(i) + Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
              end do

           end do
        end do
     end do
  end do
  do p=nC+1,nBas-nR
     do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
           do i=nC+1,nO
              
              eps = e(p) + e(i) - e(a) - e(b)
              num = ERI(p,i,a,b)**2
              
              Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
              Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
              
           end do
           do i=nC+1,nO
              
              do m=1,nS
                 num = ERI(p,a,b,i) * rho(a,i,m) * rho(b,p,m)
                 dem1 = e(p) + e(i) - e(a) - e(b)
                 dem2 = e(p) - e(b) - Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2 - num/dem1/dem2/dem2

                 num = - ERI(p,a,b,i) * rho(a,i,m) * rho(b,p,m)
                 dem1 = e(p) + e(i) - e(a) - e(b)
                 dem2 = e(a) - e(i) + Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
                 
                 num = - ERI(p,i,b,a) * rho(i,a,m) * rho(b,p,m)
                 dem1 = e(p) - e(b) - Om(m)
                 dem2 = e(a) - e(i) + Om(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
              end do

           end do
        end do
     end do
  end do

  
! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      do m=1,nS

        eps = e(a) - e(i) + Om(m)
        num = rho(a,i,m)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
