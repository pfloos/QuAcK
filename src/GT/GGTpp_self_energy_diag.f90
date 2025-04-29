subroutine GGTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOO,nVV,e,Om1,rho1,Om2,rho2,EcGM,Sig,Z,ERI)

! Compute diagonal of the correlation part of the T-matrix self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: i,j,k,a,b,c,p,m,cd,kl
  double precision              :: num,eps,dem1,dem2

! Output variables

  double precision,intent(inout)  :: EcGM
  double precision,intent(inout)  :: Sig(nBas)
  double precision,intent(inout)  :: Z(nBas)

! Initialization

  Sig(:) = 0d0
  Z(:)   = 0d0
  EcGM   = 0d0

!--------------------------------------!
! Occupied part of the Tpp self-energy !
!--------------------------------------!

!   do p=nC+1,nBas-nR
!     do i=nC+1,nO

!       do cd=1,nVV
!         eps = e(p) + e(i) - Om1(cd)
!         num = rho1(p,i,cd)**2
!         Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
!         Z(p)   = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
!       end do

!     end do
!   end do

! !------------------------------------------!
! ! Virtual part of the T-matrix self-energy !
! !------------------------------------------!

!   do p=nC+1,nBas-nR
!     do a=nO+1,nBas-nR

!       do kl=1,nOO
!         eps = e(p) + e(a) - Om2(kl)
!         num = rho2(p,a,kl)**2
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
              num = 0.5d0*(ERI(p,a,i,j) - ERI(p,a,j,i))**2

              Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
              Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2

           end do
           do a=nO+1,nBas-nR

              do m=1,nVV
                 num = - ERI(p,a,i,j) * rho1(p,a,m) * rho1(i,j,m)
                 dem1 = e(p) + e(a) - e(i) - e(j)
                 dem2 = Om1(m) - e(i) - e(j)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
              end do
              
              do m=1,nOO
                 num = - ERI(p,a,i,j) * rho2(p,a,m) * rho2(i,j,m)
                 dem1 = e(p) + e(a) - e(i) - e(j)
                 dem2 = e(p) + e(a) - Om2(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2 - num/dem1/dem2/dem2
              end do

           end do
           ! do k=nC+1,nO

           !    do m=1,nVV
           !       num = - ERI(p,i,j,k) * rho1(p,i,m) * rho1(j,k,m)
           !       dem1 = e(p) + e(i) - Om1(m)
           !       dem2 = Om1(m) - e(j) - e(k)
           !       Sig(p) = Sig(p) + num/dem1/dem2
           !       Z(p)   = Z(p)   - num/dem1/dem1/dem2
           !    end do

           ! end do
        end do
     end do
  end do
  do p=nC+1,nBas-nR
     do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR
           do i=nC+1,nO
              
              eps = e(p) + e(i) - e(a) - e(b)
              num = 0.5d0*(ERI(p,i,a,b) - ERI(p,i,b,a))**2
              
              Sig(p) = Sig(p) + num*eps/(eps**2 + eta**2)
              Z(p)    = Z(p)   - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
              
           end do
           do i=nC+1,nO
              
              do m=1,nVV
                 num = ERI(p,i,a,b) * rho1(p,i,m) * rho1(a,b,m)
                 dem1 = e(p) + e(i) - e(a) - e(b)
                 dem2 = e(p) + e(i) - Om1(m)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2 - num/dem1/dem2/dem2
              end do
              
              do m=1,nOO
                 num = ERI(p,i,a,b) * rho2(p,i,m) * rho2(a,b,m)
                 dem1 = e(p) + e(i) - e(a) - e(b)
                 dem2 = Om2(m) - e(a) - e(b)
                 Sig(p) = Sig(p) + num/dem1/dem2
                 Z(p)   = Z(p)   - num/dem1/dem1/dem2
              end do

           end do
           ! do c=nO+1,nBas-nR
           !    do m=1,nOO
           !       num = ERI(p,a,b,c) * rho2(p,a,m) * rho2(b,c,m)
           !       dem1 = e(p) + e(a) - Om2(m)
           !       dem2 = Om2(m) - e(b) - e(c)
           !       Sig(p) = Sig(p) + num/dem1/dem2
           !       Z(p)   = Z(p)   - num/dem1/dem1/dem2
           !    end do

           ! end do
        end do
     end do
  end do
  
!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  do i=nC+1,nO
    do j=nC+1,nO

      do cd=1,nVV
        eps = e(i) + e(j) - Om1(cd)
        num = 0.5d0*rho1(i,j,cd)**2
        EcGM = EcGM + num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  do a=nO+1,nBas-nR
    do b=nO+1,nBas-nR

      do kl=1,nOO
        eps = e(a) + e(b) - Om2(kl)
        num = 0.5d0*rho2(a,b,kl)**2
        EcGM = EcGM - num*eps/(eps**2 + eta**2)
      end do

    end do
  end do

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
