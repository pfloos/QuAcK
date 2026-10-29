subroutine R_Parquet_self_energy_alt_diag(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
                                          eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,          &
                                          ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,          &
                                          hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,          &
                                          EcGM,SigC,Z)

! Denominator-safe spin-adapted parquet self-energy.
! Same call signature as R_Parquet_self_energy_diag.
!
! The same-sector denominators Delta-Omega and
! Omega-e_i-e_j / Omega-e_a-e_b have been analytically removed.
!
! Exact equivalence with the original expression is obtained in the
! unregularized limit.  At finite eta this corresponds to applying
! the SRG regulator after the algebraic rearrangement.

  implicit none
  include 'parameters.h'

  double precision,intent(in)  :: eta
  integer,intent(in)           :: nOrb,nC,nO,nV,nR,nS
  integer,intent(in)           :: nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)  :: eQP(nOrb)
  double precision,intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)  :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)  :: eh_sing_Om(nS)
  double precision,intent(in)  :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)  :: eh_trip_Om(nS)

  double precision,intent(in)  :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)  :: ee_sing_Om(nVVs)
  double precision,intent(in)  :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)  :: ee_trip_Om(nVVt)

  double precision,intent(in)  :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)  :: hh_sing_Om(nOOs)
  double precision,intent(in)  :: hh_trip_rho(nOrb,nOrb,nOOt)
  double precision,intent(in)  :: hh_trip_Om(nOOt)

  double precision,intent(out) :: EcGM
  double precision,intent(out) :: SigC(nOrb)
  double precision,intent(out) :: Z(nOrb)

  call R_Parquet_self_energy_alt_omega(.true.,.false.,0,eta,eQP, &
         nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI, &
         eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om, &
         ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om, &
         hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,EcGM,SigC,Z)

end subroutine R_Parquet_self_energy_alt_diag


subroutine R_Parquet_self_energy_alt_omega(reg_PA,exact_derivative,p_eval,eta,omega, &
                                          nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,eQP,ERI,      &
                                          eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,          &
                                          ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,          &
                                          hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,          &
                                          EcGM,SigC,Z)

! Denominator-safe spin-adapted parquet self-energy.
! External frequencies omega are independent of the frozen internal eQP.
!
! The same-sector denominators Delta-Omega and
! Omega-e_i-e_j / Omega-e_a-e_b have been analytically removed.
!
! Exact equivalence with the original expression is obtained in the
! unregularized limit.  At finite eta this corresponds to applying
! the SRG regulator after the algebraic rearrangement.

  implicit none
  include 'parameters.h'

! omega varies only the external frequency; eQP and all poles stay frozen.
! p_eval=0 evaluates all active orbitals; otherwise only p_eval is evaluated.
! exact_derivative differentiates the SRG factor as well as its denominator.
! The legacy wrapper retains its previous derivative convention.
  logical,intent(in)           :: reg_PA,exact_derivative
  integer,intent(in)           :: p_eval
  double precision,intent(in)  :: eta
  integer,intent(in)           :: nOrb,nC,nO,nV,nR,nS
  integer,intent(in)           :: nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)  :: omega(nOrb),eQP(nOrb)
  double precision,intent(in)  :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)  :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)  :: eh_sing_Om(nS)
  double precision,intent(in)  :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)  :: eh_trip_Om(nS)

  double precision,intent(in)  :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)  :: ee_sing_Om(nVVs)
  double precision,intent(in)  :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)  :: ee_trip_Om(nVVt)

  double precision,intent(in)  :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)  :: hh_sing_Om(nOOs)
  double precision,intent(in)  :: hh_trip_rho(nOrb,nOrb,nOOt)
  double precision,intent(in)  :: hh_trip_Om(nOOt)

  double precision,intent(out) :: EcGM
  double precision,intent(out) :: SigC(nOrb)
  double precision,intent(out) :: Z(nOrb)

  integer                      :: p,i,j,k,a,b,c,n,p_first,p_last
  double precision             :: eps,d1,d2,num
  double precision             :: r1,r2,s1,s2
  double precision             :: dSig(nOrb)

  p_first=nC+1
  p_last=nOrb-nR
  if(p_eval /= 0) then
    if(p_eval < p_first .or. p_eval > p_last) error stop 'Invalid parquet evaluation orbital'
    p_first=p_eval
    p_last=p_eval
  end if
  SigC(:) = 0d0
  dSig(:) = 0d0
  EcGM    = 0d0

! ============================================================
! GF2
! ============================================================

  do p=p_first,p_last

    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          eps = omega(p) + eQP(a) - eQP(i) - eQP(j)

          num = 2d0*ERI(p,a,j,i)*ERI(j,i,p,a) &
              -      ERI(p,a,j,i)*ERI(j,i,a,p)

          call reg_factors(eps,r1,s1)

          SigC(p) = SigC(p) + num*r1
          dSig(p) = dSig(p) - num*s1

        end do
      end do
    end do

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          eps = omega(p) + eQP(i) - eQP(a) - eQP(b)

          num = 2d0*ERI(p,i,b,a)*ERI(b,a,p,i) &
              -      ERI(p,i,b,a)*ERI(b,a,i,p)

          call reg_factors(eps,r1,s1)

          SigC(p) = SigC(p) + num*r1
          dSig(p) = dSig(p) - num*s1

        end do
      end do
    end do

  end do

! ============================================================
! Singlet eh
! ============================================================

  do p=p_first,p_last
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do n=1,nS

          do j=nC+1,nO

            num = (0.5d0*ERI(p,a,i,j)-ERI(p,a,j,i)) &
                * eh_sing_rho(i,a,n)*eh_sing_rho(p,j,n)

            d1 = omega(p)-eQP(j)+eh_sing_Om(n)
            d2 = omega(p)-eQP(i)-eQP(j)+eQP(a)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

            num = (0.5d0*ERI(p,i,a,j)-ERI(p,i,j,a)) &
                * eh_sing_rho(a,i,n)*eh_sing_rho(p,j,n)

            d1 = eQP(a)-eQP(i)+eh_sing_Om(n)
            d2 = omega(p)-eQP(j)+eh_sing_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = (0.5d0*ERI(p,a,i,j)-ERI(p,a,j,i)) &
                * eh_sing_rho(a,i,n)*eh_sing_rho(j,p,n)

            d1 = eQP(a)-eQP(i)+eh_sing_Om(n)
            d2 = omega(p)+eQP(a)-eQP(i)-eQP(j)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

          do b=nO+1,nOrb-nR

            num = -(0.5d0*ERI(p,i,a,b)-ERI(p,i,b,a)) &
                * eh_sing_rho(i,a,n)*eh_sing_rho(b,p,n)

            d1 = omega(p)-eQP(b)-eh_sing_Om(n)
            d2 = omega(p)+eQP(i)-eQP(a)-eQP(b)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

            num = (0.5d0*ERI(p,a,i,b)-ERI(p,a,b,i)) &
                * eh_sing_rho(a,i,n)*eh_sing_rho(b,p,n)

            d1 = eQP(a)-eQP(i)+eh_sing_Om(n)
            d2 = omega(p)-eQP(b)-eh_sing_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = (0.5d0*ERI(p,i,a,b)-ERI(p,i,b,a)) &
                * eh_sing_rho(a,i,n)*eh_sing_rho(p,b,n)

            d1 = eQP(a)-eQP(i)+eh_sing_Om(n)
            d2 = omega(p)+eQP(i)-eQP(a)-eQP(b)

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
! Triplet eh
! ============================================================

  do p=p_first,p_last
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do n=1,nS

          do j=nC+1,nO

            num = 1.5d0*ERI(p,a,i,j) &
                * eh_trip_rho(i,a,n)*eh_trip_rho(p,j,n)

            d1 = omega(p)-eQP(j)+eh_trip_Om(n)
            d2 = omega(p)-eQP(i)-eQP(j)+eQP(a)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

            num = 1.5d0*ERI(p,i,a,j) &
                * eh_trip_rho(a,i,n)*eh_trip_rho(p,j,n)

            d1 = eQP(a)-eQP(i)+eh_trip_Om(n)
            d2 = omega(p)-eQP(j)+eh_trip_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = 1.5d0*ERI(p,a,i,j) &
                * eh_trip_rho(a,i,n)*eh_trip_rho(j,p,n)

            d1 = eQP(a)-eQP(i)+eh_trip_Om(n)
            d2 = omega(p)+eQP(a)-eQP(i)-eQP(j)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

          end do

          do b=nO+1,nOrb-nR

            num = -1.5d0*ERI(p,i,a,b) &
                * eh_trip_rho(i,a,n)*eh_trip_rho(b,p,n)

            d1 = omega(p)-eQP(b)-eh_trip_Om(n)
            d2 = omega(p)+eQP(i)-eQP(a)-eQP(b)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*(s1*r2+r1*s2)

            num = 1.5d0*ERI(p,a,i,b) &
                * eh_trip_rho(a,i,n)*eh_trip_rho(b,p,n)

            d1 = eQP(a)-eQP(i)+eh_trip_Om(n)
            d2 = omega(p)-eQP(b)-eh_trip_Om(n)

            call reg_factors(d1,r1,s1)
            call reg_factors(d2,r2,s2)

            SigC(p) = SigC(p) + num*r1*r2
            dSig(p) = dSig(p) - num*r1*s2

            num = 1.5d0*ERI(p,i,a,b) &
                * eh_trip_rho(a,i,n)*eh_trip_rho(p,b,n)

            d1 = eQP(a)-eQP(i)+eh_trip_Om(n)
            d2 = omega(p)+eQP(i)-eQP(a)-eQP(b)

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
! Singlet pp
! ============================================================

  call add_pp_channel(0.5d0,nVVs,nOOs,ee_sing_rho,ee_sing_Om, &
                      hh_sing_rho,hh_sing_Om)

! ============================================================
! Triplet pp
! ============================================================

  call add_pp_channel(1.5d0,nVVt,nOOt,ee_trip_rho,ee_trip_Om, &
                      hh_trip_rho,hh_trip_Om)

  Z(:) = 1d0/(1d0-dSig(:))

contains

  subroutine reg_factors(x,r,s)

    implicit none

    double precision,intent(in)  :: x
    double precision,intent(out) :: r,s
    double precision             :: reg,t,decay

    if(.not.reg_PA) then
      if(abs(x) < tiny(1d0)) error stop 'Unregularized parquet denominator is zero'
      r=1d0/x
      s=1d0/(x*x)
    else
      t=2d0*eta*x*x
      if(t < 1d-5) then
        reg=t*(1d0-t/2d0+t*t/6d0-t*t*t/24d0)
        r=2d0*eta*x*(1d0-t/2d0+t*t/6d0-t*t*t/24d0)
        s=2d0*eta*(1d0-t/2d0+t*t/6d0-t*t*t/24d0)
      else
        reg=1d0-exp(-t)
        r=reg/x
        s=reg/(x*x)
      end if
      if(exact_derivative) then
        decay=exp(-t)
        s=s-4d0*eta*decay
      end if
    end if

  end subroutine reg_factors


  subroutine add_pp_channel(weight,nVVx,nOOx,ee_rho_x,ee_Om_x,hh_rho_x,hh_Om_x)

    implicit none

    double precision,intent(in) :: weight
    integer,intent(in)          :: nVVx,nOOx
    double precision,intent(in) :: ee_rho_x(nOrb,nOrb,nVVx)
    double precision,intent(in) :: ee_Om_x(nVVx)
    double precision,intent(in) :: hh_rho_x(nOrb,nOrb,nOOx)
    double precision,intent(in) :: hh_Om_x(nOOx)

    integer                     :: pp,ii,jj,kk,aa,bb,cc,nn
    double precision            :: dd1,dd2,nnn
    double precision            :: rr1,rr2,ss1,ss2

    do pp=p_first,p_last

!     ee modes: opposite-sector pieces

      do ii=nC+1,nO
        do jj=nC+1,nO
          do nn=1,nVVx

            do kk=nC+1,nO

              nnn = -weight*ERI(pp,kk,ii,jj) &
                  * ee_rho_x(ii,jj,nn)*ee_rho_x(pp,kk,nn)

              dd1 = ee_Om_x(nn)-eQP(ii)-eQP(jj)
              dd2 = omega(pp)+eQP(kk)-ee_Om_x(nn)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*rr1*ss2

            end do

            do cc=nO+1,nOrb-nR

              nnn = -weight*ERI(pp,cc,ii,jj) &
                  * ee_rho_x(ii,jj,nn)*ee_rho_x(pp,cc,nn)

              dd1 = ee_Om_x(nn)-eQP(ii)-eQP(jj)
              dd2 = omega(pp)+eQP(cc)-eQP(ii)-eQP(jj)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*rr1*ss2

            end do

          end do
        end do
      end do

!     hh same-sector contribution: auxiliary denominator eliminated

      do ii=nC+1,nO
        do jj=nC+1,nO
          do nn=1,nOOx
            do aa=nO+1,nOrb-nR

              nnn = -weight*ERI(pp,aa,ii,jj) &
                  * hh_rho_x(ii,jj,nn)*hh_rho_x(pp,aa,nn)

              dd1 = omega(pp)+eQP(aa)-hh_Om_x(nn)
              dd2 = omega(pp)+eQP(aa)-eQP(ii)-eQP(jj)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*(ss1*rr2+rr1*ss2)

            end do
          end do
        end do
      end do

!     hh opposite-sector pieces

      do aa=nO+1,nOrb-nR
        do bb=nO+1,nOrb-nR
          do nn=1,nOOx

            do cc=nO+1,nOrb-nR

              nnn = weight*ERI(pp,cc,aa,bb) &
                  * hh_rho_x(aa,bb,nn)*hh_rho_x(pp,cc,nn)

              dd1 = hh_Om_x(nn)-eQP(aa)-eQP(bb)
              dd2 = omega(pp)+eQP(cc)-hh_Om_x(nn)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*rr1*ss2

            end do

            do kk=nC+1,nO

              nnn = weight*ERI(pp,kk,aa,bb) &
                  * hh_rho_x(aa,bb,nn)*hh_rho_x(pp,kk,nn)

              dd1 = hh_Om_x(nn)-eQP(aa)-eQP(bb)
              dd2 = omega(pp)+eQP(kk)-eQP(aa)-eQP(bb)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*rr1*ss2

            end do

          end do
        end do
      end do

!     ee same-sector addition contribution

      do aa=nO+1,nOrb-nR
        do bb=nO+1,nOrb-nR
          do nn=1,nVVx
            do ii=nC+1,nO

              nnn = weight*ERI(pp,ii,aa,bb) &
                  * ee_rho_x(aa,bb,nn)*ee_rho_x(pp,ii,nn)

              dd1 = omega(pp)+eQP(ii)-eQP(aa)-eQP(bb)
              dd2 = omega(pp)+eQP(ii)-ee_Om_x(nn)

              call reg_factors(dd1,rr1,ss1)
              call reg_factors(dd2,rr2,ss2)

              SigC(pp) = SigC(pp) + nnn*rr1*rr2
              dSig(pp) = dSig(pp) - nnn*(ss1*rr2+rr1*ss2)

            end do
          end do
        end do
      end do

    end do

  end subroutine add_pp_channel

end subroutine R_Parquet_self_energy_alt_omega
