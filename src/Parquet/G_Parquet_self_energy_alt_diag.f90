subroutine G_Parquet_self_energy_alt_diag(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eQP,ERI, &
                                          eh_rho,eh_Om,ee_rho,ee_Om,hh_rho,hh_Om,   &
                                          EcGM,SigC,Z)

! Denominator-safe form of the generalized parquet self-energy.
! Same calling sequence as G_Parquet_self_energy_diag.
!
! This is algebraically equivalent to the unregularized parquet
! self-energy and has the same safe denominator structure as the
! present generalized implementation.

  implicit none
  include 'parameters.h'

! Input

  double precision,intent(in)  :: eta
  integer,intent(in)           :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)  :: eQP(nOrb)
  double precision,intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)  :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)  :: eh_Om(nS)
  double precision,intent(in)  :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)  :: ee_Om(nVV)
  double precision,intent(in)  :: hh_rho(nOrb,nOrb,nOO)
  double precision,intent(in)  :: hh_Om(nOO)

! Output

  double precision,intent(out) :: EcGM
  double precision,intent(out) :: SigC(nOrb)
  double precision,intent(out) :: Z(nOrb)

! Local

  integer                      :: p,i,j,k,a,b,c,n
  double precision             :: eps,d1,d2,num
  double precision             :: r1,r2,s1,s2
  double precision             :: dSig(nOrb)

  SigC(:) = 0d0
  dSig(:) = 0d0
  EcGM    = 0d0

! ============================================================
! GF2
! ============================================================

  do p=nC+1,nOrb-nR

    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          eps = eQP(p) + eQP(a) - eQP(i) - eQP(j)
          num = 0.5d0*(ERI(p,a,j,i)-ERI(p,a,i,j))**2

          call reg_factors(eps,r1,s1)
          SigC(p) = SigC(p) + num*r1
          dSig(p) = dSig(p) - num*s1

        end do
      end do
    end do

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          eps = eQP(p) + eQP(i) - eQP(a) - eQP(b)
          num = 0.5d0*(ERI(p,i,b,a)-ERI(p,i,a,b))**2

          call reg_factors(eps,r1,s1)
          SigC(p) = SigC(p) + num*r1
          dSig(p) = dSig(p) - num*s1

        end do
      end do
    end do

  end do

! ============================================================
! eh channel
! ============================================================

  do p=nC+1,nOrb-nR
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do n=1,nS

!       Removal branch

          do j=nC+1,nO

!           Same-sector contribution:
!           only genuine bare and dressed one-body poles.

            num = (ERI(p,a,i,j)-ERI(p,a,j,i)) &
                * eh_rho(i,a,n)*eh_rho(p,j,n)

            d1 = eQP(p) - eQP(j) + eh_Om(n)
            d2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

!           Backward-BSE component; Delta + Omega is harmless
!           for a stable Aufbau reference.

            num = (ERI(p,i,a,j)-ERI(p,i,j,a)) &
                * eh_rho(a,i,n)*eh_rho(p,j,n)

            d1 = eQP(a) - eQP(i) + eh_Om(n)
            d2 = eQP(p) - eQP(j) + eh_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = (ERI(p,a,i,j)-ERI(p,a,j,i)) &
                * eh_rho(a,i,n)*eh_rho(j,p,n)

            d1 = eQP(a) - eQP(i) + eh_Om(n)
            d2 = eQP(p) - eQP(i) - eQP(j) + eQP(a)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

!       Addition branch

          do b=nO+1,nOrb-nR

            num = -(ERI(p,i,a,b)-ERI(p,i,b,a)) &
                * eh_rho(i,a,n)*eh_rho(b,p,n)

            d1 = eQP(p) - eQP(b) - eh_Om(n)
            d2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

            num = (ERI(p,a,i,b)-ERI(p,a,b,i)) &
                * eh_rho(a,i,n)*eh_rho(b,p,n)

            d1 = eQP(a) - eQP(i) + eh_Om(n)
            d2 = eQP(p) - eQP(b) - eh_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = (ERI(p,i,a,b)-ERI(p,i,b,a)) &
                * eh_rho(a,i,n)*eh_rho(p,b,n)

            d1 = eQP(a) - eQP(i) + eh_Om(n)
            d2 = eQP(p) + eQP(i) - eQP(a) - eQP(b)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

        end do
      end do
    end do
  end do

! ============================================================
! pp channel
! ============================================================

  do p=nC+1,nOrb-nR

!   Removal-like terms

    do i=nC+1,nO
      do j=nC+1,nO

        do n=1,nVV

          do k=nC+1,nO

            num = -ERI(p,k,i,j)*ee_rho(i,j,n)*ee_rho(p,k,n)

            d1 = ee_Om(n) - eQP(i) - eQP(j)
            d2 = eQP(p) + eQP(k) - ee_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

          do c=nO+1,nOrb-nR

            num = -ERI(p,c,i,j)*ee_rho(i,j,n)*ee_rho(p,c,n)

            d1 = ee_Om(n) - eQP(i) - eQP(j)
            d2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

        end do

        do n=1,nOO
          do c=nO+1,nOrb-nR

            num = -ERI(p,c,i,j)*hh_rho(i,j,n)*hh_rho(p,c,n)

            d1 = eQP(p) + eQP(c) - hh_Om(n)
            d2 = eQP(p) + eQP(c) - eQP(i) - eQP(j)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

          end do
        end do

      end do
    end do

!   Addition-like terms

    do a=nO+1,nOrb-nR
      do b=nO+1,nOrb-nR

        do n=1,nOO

          do c=nO+1,nOrb-nR

            num = ERI(p,c,a,b)*hh_rho(a,b,n)*hh_rho(p,c,n)

            d1 = hh_Om(n) - eQP(a) - eQP(b)
            d2 = eQP(p) + eQP(c) - hh_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

          do k=nC+1,nO

            num = ERI(p,k,a,b)*hh_rho(a,b,n)*hh_rho(p,k,n)

            d1 = hh_Om(n) - eQP(a) - eQP(b)
            d2 = eQP(p) + eQP(k) - eQP(a) - eQP(b)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

        end do

        do n=1,nVV
          do k=nC+1,nO

            num = ERI(p,k,a,b)*ee_rho(a,b,n)*ee_rho(p,k,n)

            d1 = eQP(p) + eQP(k) - eQP(a) - eQP(b)
            d2 = eQP(p) + eQP(k) - ee_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

          end do
        end do

      end do
    end do

  end do

  Z(:) = 1d0/(1d0-dSig(:))

contains

  subroutine reg_factors(x,r,s)

    implicit none

    double precision,intent(in)  :: x
    double precision,intent(out) :: r,s
    double precision             :: reg

    reg = 1d0 - exp(-2d0*eta*x*x)

    r = reg/x
    s = reg/(x*x)

  end subroutine reg_factors

end subroutine G_Parquet_self_energy_alt_diag
