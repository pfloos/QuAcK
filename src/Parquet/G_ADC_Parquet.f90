subroutine G_ADC_Parquet(reg_PA,eta_1b,                              &
                                nOrb,nC,nO,nV,nR,nS,nOO,nVV,                 &
                                eHF,eQP,ERI,                                  &
                                eh_Om,eh_rho,                                 &
                                ee_Om,ee_rho,hh_Om,hh_rho,                    &
                                eADC,ZADC)

! Full Hermitian / PSD ADC-Parquet construction in spin-orbital form.
!
! The ADC gauge is controlled by six relative weights defined in the
! "ADC gauge controls" block below.  The code normalizes these weights so
! that
!
!   alpha_B-**2 + 2 alpha_eh-**2 + alpha_hh-**2 = 1
!   alpha_B+**2 + 2 alpha_eh+**2 + alpha_ee+**2 = 1
!
! and then enforces the parent-self-energy matching conditions
!
!   beta_gamma = 1/(2 alpha_gamma)
!   C_gamma,delta = Xi_gamma,delta/(2 alpha_gamma alpha_delta)
!
! Thus changing the gauge only requires changing the six weights below.
! The default choice g=1 for every sector gives the present symmetric gauge
!
!   alpha = 1/2,  beta = 1,  C = 2 Xi.
!
! The native parquet self-energy is distributed as follows:
!
!   * GF2 + one-bare-pole pieces          -> U_B
!   * one-eh-pole pieces                  -> U_eh
!   * one-hh/ee-pole pieces               -> U_hh / U_ee
!   * B <-> eh mixed-pole pieces          -> C_B_eh
!   * B <-> hh / ee mixed-pole pieces     -> C_B_hh / C_B_ee
!
! Internal denominators entering the D corrections can be regularized
! with the same one-body SRG prescription as the parent parquet
! self-energy:
!
!   1/Delta -> [1 - exp(-2 eta_1b Delta**2)]/Delta.
!
! The regularization is controlled by reg_PA and eta_1b.  The genuine
! ADC pole energies in the K blocks and the denominator-free C blocks
! are left unchanged.

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: reg_PA
  double precision,intent(in)   :: eta_1b
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

! ADC gauge variables

  double precision              :: gBm,gEhm,gHhm,gBp,gEhp,gEep
  double precision              :: normMinus,normPlus
  double precision              :: alphaBm,alphaEhm,alphaHhm
  double precision              :: alphaBp,alphaEhp,alphaEep
  double precision              :: betaBm,betaEhm,betaHhm
  double precision              :: betaBp,betaEhp,betaEep
  double precision              :: cBmEhm,cBmHhm,cBpEhp,cBpEep
  double precision,parameter    :: gauge_tol = 1d-12

  nOcc = nO - nC
  nVir = nOrb - nR - nO
  nAct = nOrb - nR - nC

! ======================================================================
! ADC gauge controls
! ======================================================================
!
! Edit ONLY these six relative weights to explore the gauge freedom.
! An overall rescaling of all weights within one branch is irrelevant,
! because the weights are normalized below.
!
! Symmetric gauge: all six weights = 1.
!
! Example: emphasize the bare B sectors relative to the dressed sectors:
!   gBm = 2d0 ; gEhm = 1d0 ; gHhm = 1d0
!   gBp = 2d0 ; gEhp = 1d0 ; gEep = 1d0
!
! None of the normalized alpha values may vanish, because the exact parent
! matching then requires beta or C to diverge.

  gBm  = 1d0
  gEhm = 1d0
  gHhm = 1d0

  gBp  = 1d0
  gEhp = 1d0
  gEep = 1d0

! Normalize the relative weights so that GF2 is reproduced exactly once.

  normMinus = sqrt(gBm*gBm + 2d0*gEhm*gEhm + gHhm*gHhm)
  normPlus  = sqrt(gBp*gBp + 2d0*gEhp*gEhp + gEep*gEep)

  if(normMinus < gauge_tol .or. normPlus < gauge_tol) then
    write(*,*) 'ERROR in G_ADC_Parquet: vanishing ADC gauge norm.'
    stop
  end if

  alphaBm  = gBm  / normMinus
  alphaEhm = gEhm / normMinus
  alphaHhm = gHhm / normMinus

  alphaBp  = gBp  / normPlus
  alphaEhp = gEhp / normPlus
  alphaEep = gEep / normPlus

  if(abs(alphaBm ) < gauge_tol .or. abs(alphaEhm) < gauge_tol .or. &
     abs(alphaHhm) < gauge_tol .or. abs(alphaBp ) < gauge_tol .or. &
     abs(alphaEhp) < gauge_tol .or. abs(alphaEep) < gauge_tol) then
    write(*,*) 'ERROR in G_ADC_Parquet: all ADC alpha values must be nonzero.'
    stop
  end if

! Parent-one-pole matching: 2 alpha_gamma beta_gamma = 1.

  betaBm  = 1d0/(2d0*alphaBm )
  betaEhm = 1d0/(2d0*alphaEhm)
  betaHhm = 1d0/(2d0*alphaHhm)

  betaBp  = 1d0/(2d0*alphaBp )
  betaEhp = 1d0/(2d0*alphaEhp)
  betaEep = 1d0/(2d0*alphaEep)

! Parent-mixed-pole matching:
! C_gamma,delta = Xi_gamma,delta/(2 alpha_gamma alpha_delta).

  cBmEhm = 1d0/(2d0*alphaBm*alphaEhm)
  cBmHhm = 1d0/(2d0*alphaBm*alphaHhm)
  cBpEhp = 1d0/(2d0*alphaBp*alphaEhp)
  cBpEep = 1d0/(2d0*alphaBp*alphaEep)

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
  write(*,'(1X,A34,L4)') 'SRG regularization enabled =',reg_PA
  if(reg_PA) write(*,'(1X,A34,1X,ES12.4)') 'One-body SRG parameter eta_1b =',eta_1b
  write(*,*)

  write(*,*)'ADC gauge parameters:'
  write(*,'(1X,A8,3(1X,F12.8))') 'alpha - ',alphaBm,alphaEhm,alphaHhm
  write(*,'(1X,A8,3(1X,F12.8))') 'beta  - ',betaBm,betaEhm,betaHhm
  write(*,'(1X,A8,3(1X,F12.8))') 'alpha + ',alphaBp,alphaEhp,alphaEep
  write(*,'(1X,A8,3(1X,F12.8))') 'beta  + ',betaBp,betaEhp,betaEep
  write(*,'(1X,A12,2(1X,F12.8))') 'C scales - ',cBmEhm,cBmHhm
  write(*,'(1X,A12,2(1X,F12.8))') 'C scales + ',cBpEhp,cBpEep
  write(*,*)

! ======================================================================
! Physical block
! ======================================================================

  do p=1,nOrb
    H(p,p) = eHF(p)
  end do

! ======================================================================
! B- block: U_B- = alpha_B- V- + beta_B- D_B-
! ======================================================================

  ib = oBm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(i) + eQP(j) - eQP(a)

        do p=nC+1,nOrb-nR
          val = alphaBm*vminus(p,a,i,j) + betaBm*DBminus(i,j,a,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! eh- block: U_eh- = alpha_eh- rho + beta_eh- D_eh-
! ======================================================================

  ib = oEhm

  do j=nC+1,nO
    do n=1,nS

      ib = ib + 1
      H(ib,ib) = eQP(j) - eh_Om(n)

      do p=nC+1,nOrb-nR
        val = alphaEhm*eh_rho(p,j,n) + betaEhm*Dehminus(j,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! hh- block: U_hh- = alpha_hh- rho + beta_hh- D_hh-
! ======================================================================

  ib = oHhm

  do a=nO+1,nOrb-nR
    do n=1,nOO

      ib = ib + 1
      H(ib,ib) = hh_Om(n) - eQP(a)

      do p=nC+1,nOrb-nR
        val = alphaHhm*hh_rho(p,a,n) + betaHhm*Dhhminus(a,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! B+ block: U_B+ = alpha_B+ V+ + beta_B+ D_B+
! ======================================================================

  ib = oBp

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(a) + eQP(b) - eQP(i)

        do p=nC+1,nOrb-nR
          val = alphaBp*vplus(p,i,a,b) + betaBp*DBplus(i,a,b,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! eh+ block: U_eh+ = alpha_eh+ rho + beta_eh+ D_eh+
! ======================================================================

  ib = oEhp

  do b=nO+1,nOrb-nR
    do n=1,nS

      ib = ib + 1
      H(ib,ib) = eQP(b) + eh_Om(n)

      do p=nC+1,nOrb-nR
        val = alphaEhp*eh_rho(b,p,n) + betaEhp*Dehplus(b,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! ee+ block: U_ee+ = alpha_ee+ rho + beta_ee+ D_ee+
! ======================================================================

  ib = oEep

  do i=nC+1,nO
    do n=1,nVV

      ib = ib + 1
      H(ib,ib) = ee_Om(n) - eQP(i)

      do p=nC+1,nOrb-nR
        val = alphaEep*ee_rho(p,i,n) + betaEep*Deeplus(i,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! C_B-,eh- = cBmEhm Xi_B-,eh-
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

            if(j == k) val = val + cBmEhm*eh_rho(i,a,n)
            if(i == k) val = val - cBmEhm*eh_rho(j,a,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B-,hh- = cBmHhm Xi_B-,hh-
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

            if(a == b) val = -cBmHhm*hh_rho(i,j,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B+,eh+ = cBpEhp Xi_B+,eh+
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

            if(a == c) val = val + cBpEhp*eh_rho(i,b,n)
            if(b == c) val = val - cBpEhp*eh_rho(i,a,n)

            H(ib,jb) = val
            H(jb,ib) = val

          end do
        end do

      end do
    end do
  end do

! ======================================================================
! C_B+,ee+ = cBpEep Xi_B+,ee+
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

            if(i == j) val = cBpEep*ee_rho(a,b,n)

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
! Internal denominator regularization
! ======================================================================
!
! Use the same one-body SRG regularization as the parent parquet
! self-energy and ADC-G3W2.  Only the explicit internal denominators in
! the D corrections pass through this helper.  K and C are untouched.

  double precision function regularized_inverse(dem)

    implicit none
    double precision,intent(in) :: dem
    double precision,external   :: SRG_reg

    if(reg_PA) then
      regularized_inverse = SRG_reg(dem,eta_1b)
    else
      regularized_inverse = 1d0/dem
    end if

  end function regularized_inverse

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
      DBminus = DBminus + eh_rho(aa,ii,nn)*eh_rho(jj,pp,nn)*regularized_inverse(dem)

      dem = eQP(aa) - eQP(jj) + eh_Om(nn)
      DBminus = DBminus - eh_rho(aa,jj,nn)*eh_rho(ii,pp,nn)*regularized_inverse(dem)

    end do

!   pp/ee one-bare-pole contribution

    do nn=1,nVV

      dem = ee_Om(nn) - eQP(ii) - eQP(jj)
      DBminus = DBminus - ee_rho(ii,jj,nn)*ee_rho(pp,aa,nn)*regularized_inverse(dem)

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
      DBplus = DBplus + eh_rho(aa,ii,nn)*eh_rho(pp,bb,nn)*regularized_inverse(dem)

      dem = eQP(bb) - eQP(ii) + eh_Om(nn)
      DBplus = DBplus - eh_rho(bb,ii,nn)*eh_rho(pp,aa,nn)*regularized_inverse(dem)

    end do

!   pp/hh one-bare-pole contribution

    do nn=1,nOO

      dem = hh_Om(nn) - eQP(aa) - eQP(bb)
      DBplus = DBplus + hh_rho(aa,bb,nn)*hh_rho(pp,ii,nn)*regularized_inverse(dem)

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
          + (ERI(pp,ii,aa,jj)-ERI(pp,ii,jj,aa))*eh_rho(aa,ii,nn)*regularized_inverse(dem)

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
          + (ERI(pp,aa,ii,bb)-ERI(pp,aa,bb,ii))*eh_rho(aa,ii,nn)*regularized_inverse(dem)

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
          + (ERI(pp,aa,bb,cc)-ERI(pp,aa,cc,bb))*hh_rho(bb,cc,nn)*regularized_inverse(dem)

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
          - (ERI(pp,ii,jj,kk)-ERI(pp,ii,kk,jj))*ee_rho(jj,kk,nn)*regularized_inverse(dem)

      end do
    end do

  end function Deeplus

end subroutine G_ADC_Parquet

