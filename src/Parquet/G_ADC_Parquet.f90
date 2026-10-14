subroutine G_ADC_Parquet(nOrb,nC,nO,nV,nR,nS,nOO,nVV,                 &
                                eHF,eQP,ERI,                                  &
                                eh_Om,eh_rho,                                 &
                                ee_Om,ee_rho,hh_Om,hh_rho,                    &
                                eADC,ZADC)

! Full Hermitian / PSD ADC-Parquet construction in spin-orbital form.
!
! Symmetric ADC gauge:
!
!   alpha_B = alpha_eh = alpha_hh/ee = 1/2
!
! The native parquet self-energy is distributed as follows:
!
!   * GF2 + one-bare-pole pieces          -> U_B
!   * one-eh-pole pieces                  -> U_eh
!   * one-hh/ee-pole pieces               -> U_hh / U_ee
!   * B <-> eh mixed-pole pieces          -> C_B_eh
!   * B <-> hh / ee mixed-pole pieces     -> C_B_hh / C_B_ee
!
! The Hermitian resolvent therefore contains every native PA term and
! generates additional PSD ADC completion terms.
!
! No SRG regulator is applied inside the ADC decomposition.  The
! identities used here correspond to the unregularized parent PA
! self-energy.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eHF(nOrb),eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)
  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Output variables

  double precision,intent(out)  :: eADC(nOrb),ZADC(nOrb)

! Local variables

  integer                       :: nOcc,nVir,nAct
  integer                       :: nBm,nEhm,nHhm,nBp,nEhp,nEep,nH
  integer                       :: oBm,oEhm,oHhm,oBp,oEhp,oEep
  integer                       :: p,i,j,k,a,b,c,n,s,idx,it
  integer                       :: ib,jb
  double precision              :: val,cutoff,tol_deg,wgrp,wbest
  double precision              :: asym,sumZ
  double precision,allocatable  :: H(:,:),ePA(:),Zstate(:)

  nOcc = nO - nC
  nVir = nOrb - nR - nO
  nAct = nOrb - nR - nC

! Auxiliary dimensions

  nBm  = nOcc*(nOcc-1)*nVir/2
  nEhm = nOcc*nS
  nHhm = nVir*nOO

  nBp  = nOcc*nVir*(nVir-1)/2
  nEhp = nVir*nS
  nEep = nOcc*nVV

! Offsets

  oBm  = nOrb
  oEhm = oBm  + nBm
  oHhm = oEhm + nEhm
  oBp  = oHhm + nHhm
  oEhp = oBp  + nBp
  oEep = oEhp + nEhp

  nH = oEep + nEep

  cutoff  = 1d-2
  tol_deg = 1d-8

  allocate(H(nH,nH),ePA(nH),Zstate(nH))

  H(:,:)    = 0d0
  ePA(:)    = 0d0
  Zstate(:) = 0d0
  eADC(:)   = eHF(:)
  ZADC(:)   = 1d0

  write(*,*)
  write(*,*)'***********************************************'
  write(*,*)'* Full generalized Hermitian ADC-Parquet     *'
  write(*,*)'***********************************************'
  write(*,'(1X,A24,I10)') 'B- dimension       =',nBm
  write(*,'(1X,A24,I10)') 'eh- dimension      =',nEhm
  write(*,'(1X,A24,I10)') 'hh- dimension      =',nHhm
  write(*,'(1X,A24,I10)') 'B+ dimension       =',nBp
  write(*,'(1X,A24,I10)') 'eh+ dimension      =',nEhp
  write(*,'(1X,A24,I10)') 'ee+ dimension      =',nEep
  write(*,'(1X,A24,I10)') 'Total dimension    =',nH
  write(*,*)

! ======================================================================
! Physical block
! ======================================================================

  do p=1,nOrb
    H(p,p) = eHF(p)
  end do

! ======================================================================
! B- block: K_B- and U_B- = 1/2 V- + D_B-
! ======================================================================

  ib = oBm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(i) + eQP(j) - eQP(a)

        do p=nC+1,nOrb-nR
          val = 0.5d0*vminus(p,a,i,j) + DBminus(i,j,a,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! eh- block: K_eh- and U_eh- = 1/2 rho + D_eh-
! ======================================================================

  ib = oEhm

  do j=nC+1,nO
    do n=1,nS

      ib = ib + 1
      H(ib,ib) = eQP(j) - eh_Om(n)

      do p=nC+1,nOrb-nR
        val = 0.5d0*eh_rho(p,j,n) + Dehminus(j,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! hh- block: K_hh- and U_hh- = 1/2 rho + D_hh-
! ======================================================================

  ib = oHhm

  do a=nO+1,nOrb-nR
    do n=1,nOO

      ib = ib + 1
      H(ib,ib) = hh_Om(n) - eQP(a)

      do p=nC+1,nOrb-nR
        val = 0.5d0*hh_rho(p,a,n) + Dhhminus(a,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! B+ block
! ======================================================================

  ib = oBp

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(a) + eQP(b) - eQP(i)

        do p=nC+1,nOrb-nR
          val = 0.5d0*vplus(p,i,a,b) + DBplus(i,a,b,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! eh+ block
! ======================================================================

  ib = oEhp

  do b=nO+1,nOrb-nR
    do n=1,nS

      ib = ib + 1
      H(ib,ib) = eQP(b) + eh_Om(n)

      do p=nC+1,nOrb-nR
        val = 0.5d0*eh_rho(b,p,n) + Dehplus(b,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! ee+ block
! ======================================================================

  ib = oEep

  do i=nC+1,nO
    do n=1,nVV

      ib = ib + 1
      H(ib,ib) = ee_Om(n) - eQP(i)

      do p=nC+1,nOrb-nR
        val = 0.5d0*ee_rho(p,i,n) + Deeplus(i,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! C_B-,eh- = 2 Xi_B-,eh-
!
! Xi_{ija,k n} = delta_jk rho_ia,n - delta_ik rho_ja,n
! ======================================================================

  ib = oBm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        jb = oEhm

        do k=nC+1,nO
          do n=1,nS

            jb = jb + 1
            val = 0d0

            if(j == k) val = val + 2d0*eh_rho(i,a,n)
            if(i == k) val = val - 2d0*eh_rho(j,a,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B-,hh- = 2 Xi_B-,hh-
!
! Xi_{ija,b m} = - delta_ab rho_ij,m
! ======================================================================

  ib = oBm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        jb = oHhm

        do b=nO+1,nOrb-nR
          do n=1,nOO

            jb = jb + 1
            val = 0d0

            if(a == b) val = -2d0*hh_rho(i,j,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B+,eh+ = 2 Xi_B+,eh+
!
! Xi_{iab,c n} = delta_ac rho_ib,n - delta_bc rho_ia,n
! ======================================================================

  ib = oBp

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        ib = ib + 1
        jb = oEhp

        do c=nO+1,nOrb-nR
          do n=1,nS

            jb = jb + 1
            val = 0d0

            if(a == c) val = val + 2d0*eh_rho(i,b,n)
            if(b == c) val = val - 2d0*eh_rho(i,a,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B+,ee+ = 2 Xi_B+,ee+
!
! Xi_{iab,j m} = delta_ij rho_ab,m
! ======================================================================

  ib = oBp

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        ib = ib + 1
        jb = oEep

        do j=nC+1,nO
          do n=1,nVV

            jb = jb + 1
            val = 0d0

            if(i == j) val = 2d0*ee_rho(a,b,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! Hermiticity check
! ======================================================================

  asym = maxval(abs(H-transpose(H)))
  write(*,'(1X,A34,1X,ES12.4)') 'Max |H-H^T| before diagonalization =',asym

  H = 0.5d0*(H + transpose(H))

! ======================================================================
! Diagonalization
! ======================================================================

  call diagonalize_matrix(nH,H,ePA)

! Total active one-particle spectral weight

  do s=1,nH
    do p=nC+1,nOrb-nR
      Zstate(s) = Zstate(s) + H(p,s)**2
    end do
  end do

  sumZ = sum(Zstate)

  write(*,'(1X,A30,F16.10)') 'Sum of active 1p weights      =',sumZ
  write(*,'(1X,A30,I16)')    'Expected active 1p dimension =',nAct
  write(*,*)

! Dominant pole associated with each orbital, grouping exact degeneracies

  do p=nC+1,nOrb-nR

    wbest = -1d0
    idx   = 1
    s     = 1

    do while(s <= nH)

      it = s
      do while(it < nH)
        if(abs(ePA(it+1)-ePA(s)) > tol_deg) exit
        it = it + 1
      end do

      wgrp = sum(H(p,s:it)**2)

      if(wgrp > wbest) then
        wbest = wgrp
        idx   = s
      end if

      s = it + 1

    end do

    eADC(p) = ePA(idx)
    ZADC(p) = wbest

  end do

! Print significant poles

  write(*,*)'-------------------------------------------------------'
  write(*,*)'      Full generalized PSD ADC-Parquet spectrum         '
  write(*,*)'-------------------------------------------------------'
  write(*,'(1X,A6,2X,A16,2X,A16)') '#','energy (eV)','1p weight'

  do s=1,nH
    if(Zstate(s) > cutoff) then
      write(*,'(1X,I6,2X,F16.8,2X,F16.8)') s,ePA(s)*HaToeV,Zstate(s)
    end if
  end do

  write(*,*)'-------------------------------------------------------'
  write(*,*)

  deallocate(H,ePA,Zstate)

contains

! ======================================================================
! Antisymmetrized bare vertices
! ======================================================================

  double precision function vminus(pp,aa,ii,jj)

    implicit none
    integer,intent(in) :: pp,aa,ii,jj

    vminus = ERI(pp,aa,ii,jj) - ERI(pp,aa,jj,ii)

  end function vminus


  double precision function vplus(pp,ii,aa,bb)

    implicit none
    integer,intent(in) :: pp,ii,aa,bb

    vplus = ERI(pp,ii,aa,bb) - ERI(pp,ii,bb,aa)

  end function vplus

! ======================================================================
! B- correction
! ======================================================================

  double precision function DBminus(ii,jj,aa,pp)

    implicit none
    integer,intent(in) :: ii,jj,aa,pp

    integer            :: nn
    double precision   :: dem

    DBminus = 0d0

!   eh one-bare-pole contribution

    do nn=1,nS

      dem = eQP(aa) - eQP(ii) + eh_Om(nn)
      DBminus = DBminus + eh_rho(aa,ii,nn)*eh_rho(jj,pp,nn)/dem

      dem = eQP(aa) - eQP(jj) + eh_Om(nn)
      DBminus = DBminus - eh_rho(aa,jj,nn)*eh_rho(ii,pp,nn)/dem

    end do

!   pp/ee one-bare-pole contribution

    do nn=1,nVV

      dem = ee_Om(nn) - eQP(ii) - eQP(jj)
      DBminus = DBminus - ee_rho(ii,jj,nn)*ee_rho(pp,aa,nn)/dem

    end do

  end function DBminus

! ======================================================================
! B+ correction
! ======================================================================

  double precision function DBplus(ii,aa,bb,pp)

    implicit none
    integer,intent(in) :: ii,aa,bb,pp

    integer            :: nn
    double precision   :: dem

    DBplus = 0d0

!   eh one-bare-pole contribution

    do nn=1,nS

      dem = eQP(aa) - eQP(ii) + eh_Om(nn)
      DBplus = DBplus + eh_rho(aa,ii,nn)*eh_rho(pp,bb,nn)/dem

      dem = eQP(bb) - eQP(ii) + eh_Om(nn)
      DBplus = DBplus - eh_rho(bb,ii,nn)*eh_rho(pp,aa,nn)/dem

    end do

!   pp/hh one-bare-pole contribution

    do nn=1,nOO

      dem = hh_Om(nn) - eQP(aa) - eQP(bb)
      DBplus = DBplus + hh_rho(aa,bb,nn)*hh_rho(pp,ii,nn)/dem

    end do

  end function DBplus

! ======================================================================
! eh- correction
! ======================================================================

  double precision function Dehminus(jj,nn,pp)

    implicit none
    integer,intent(in) :: jj,nn,pp

    integer            :: ii,aa
    double precision   :: dem

    Dehminus = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR

        dem = eQP(aa) - eQP(ii) + eh_Om(nn)

        Dehminus = Dehminus                                          &
          + (ERI(pp,ii,aa,jj)-ERI(pp,ii,jj,aa))*eh_rho(aa,ii,nn)/dem

      end do
    end do

  end function Dehminus

! ======================================================================
! eh+ correction
! ======================================================================

  double precision function Dehplus(bb,nn,pp)

    implicit none
    integer,intent(in) :: bb,nn,pp

    integer            :: ii,aa
    double precision   :: dem

    Dehplus = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR

        dem = eQP(aa) - eQP(ii) + eh_Om(nn)

        Dehplus = Dehplus                                            &
          + (ERI(pp,aa,ii,bb)-ERI(pp,aa,bb,ii))*eh_rho(aa,ii,nn)/dem

      end do
    end do

  end function Dehplus

! ======================================================================
! hh- correction
! ======================================================================

  double precision function Dhhminus(aa,nn,pp)

    implicit none
    integer,intent(in) :: aa,nn,pp

    integer            :: bb,cc
    double precision   :: dem

    Dhhminus = 0d0

    do bb=nO+1,nOrb-nR
      do cc=bb+1,nOrb-nR

        dem = hh_Om(nn) - eQP(bb) - eQP(cc)

        Dhhminus = Dhhminus                                          &
          + (ERI(pp,aa,bb,cc)-ERI(pp,aa,cc,bb))*hh_rho(bb,cc,nn)/dem

      end do
    end do

  end function Dhhminus

! ======================================================================
! ee+ correction
! ======================================================================

  double precision function Deeplus(ii,nn,pp)

    implicit none
    integer,intent(in) :: ii,nn,pp

    integer            :: jj,kk
    double precision   :: dem

    Deeplus = 0d0

    do jj=nC+1,nO
      do kk=jj+1,nO

        dem = ee_Om(nn) - eQP(jj) - eQP(kk)

        Deeplus = Deeplus                                            &
          - (ERI(pp,ii,jj,kk)-ERI(pp,ii,kk,jj))*ee_rho(jj,kk,nn)/dem

      end do
    end do

  end function Deeplus

end subroutine G_ADC_Parquet
