subroutine R_ADC_Parquet(reg_PA,eta_1b,                                      &
                         nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,           &
                         eHF,eQP,ERI,                                         &
                         eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho,      &
                         ee_sing_Om,ee_sing_rho,ee_trip_Om,ee_trip_rho,      &
                         hh_sing_Om,hh_sing_rho,hh_trip_Om,hh_trip_rho,      &
                         eADC,ZADC)

! Full Hermitian / PSD ADC-Parquet construction in restricted spin-adapted form.
!
! The bare 2h1p/2p1h space is represented by the three spin-adapted doublet
! sectors used in restricted IP/EA-ADC.  The dressed eh and pp auxiliary
! sectors are resolved into singlet and triplet components.
!
! The external couplings correspond to the equal-weight ADC construction.
! Internal response denominators entering the higher-order U corrections can
! be regularized with the same one-body SRG prescription as the parent parquet
! self-energy,
!
!   1/Delta -> [1 - exp(-2 eta_1b Delta**2)]/Delta.
!
! The genuine ADC pole energies in the K blocks and the denominator-free C
! blocks are left unchanged.

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: reg_PA
  double precision,intent(in)   :: eta_1b
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  integer,intent(in)            :: nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: eHF(nOrb),eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: eh_sing_Om(nS),eh_trip_Om(nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)

  double precision,intent(in)   :: ee_sing_Om(nVVs),ee_trip_Om(nVVt)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)

  double precision,intent(in)   :: hh_sing_Om(nOOs),hh_trip_Om(nOOt)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOOt)

! Output variables

  double precision,intent(out)  :: eADC(nOrb),ZADC(nOrb)

! Local variables

  integer                       :: nOcc,nVir,nAct
  integer                       :: nIm,nIIm,nIIIm,nBm
  integer                       :: nIp,nIIp,nIIIp,nBp
  integer                       :: nEhmS,nEhmT,nHhmS,nHhmT
  integer                       :: nEhpS,nEhpT,nEepS,nEepT
  integer                       :: oBm,oEhmS,oEhmT,oHhmS,oHhmT
  integer                       :: oBp,oEhpS,oEhpT,oEepS,oEepT,nH

  integer                       :: p,i,j,k,a,b,c,n,s,idx,it
  integer                       :: ib,jb
  double precision              :: val,cutoff,tol_deg,wgrp,wbest
  double precision              :: asym,sumZ
  double precision,allocatable  :: H(:,:),ePA(:),Zstate(:)

  nOcc = nO - nC
  nVir = nOrb - nR - nO
  nAct = nOrb - nR - nC

! Bare removal sectors

  nIm   = nOcc*nVir
  nIIm  = nOcc*(nOcc-1)*nVir/2
  nIIIm = nIIm
  nBm   = nIm + nIIm + nIIIm

! Bare addition sectors

  nIp   = nOcc*nVir
  nIIp  = nOcc*nVir*(nVir-1)/2
  nIIIp = nIIp
  nBp   = nIp + nIIp + nIIIp

! Dressed sectors

  nEhmS = nOcc*nS
  nEhmT = nOcc*nS
  nHhmS = nVir*nOOs
  nHhmT = nVir*nOOt

  nEhpS = nVir*nS
  nEhpT = nVir*nS
  nEepS = nOcc*nVVs
  nEepT = nOcc*nVVt

! Offsets

  oBm   = nOrb
  oEhmS = oBm   + nBm
  oEhmT = oEhmS + nEhmS
  oHhmS = oEhmT + nEhmT
  oHhmT = oHhmS + nHhmS

  oBp   = oHhmT + nHhmT
  oEhpS = oBp   + nBp
  oEhpT = oEhpS + nEhpS
  oEepS = oEhpT + nEhpT
  oEepT = oEepS + nEepS

  nH = oEepT + nEepT

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
  write(*,*)'* Full restricted Hermitian ADC-Parquet      *'
  write(*,*)'***********************************************'
  write(*,'(1X,A34,L4)')            'SRG regularization enabled =',reg_PA
  if(reg_PA) write(*,'(1X,A34,1X,ES12.4)') 'One-body SRG parameter eta_1b =',eta_1b
  write(*,'(1X,A24,I10)') 'B- dimension       =',nBm
  write(*,'(1X,A24,I10)') 'eh- singlet        =',nEhmS
  write(*,'(1X,A24,I10)') 'eh- triplet        =',nEhmT
  write(*,'(1X,A24,I10)') 'pp- singlet        =',nHhmS
  write(*,'(1X,A24,I10)') 'pp- triplet        =',nHhmT
  write(*,'(1X,A24,I10)') 'B+ dimension       =',nBp
  write(*,'(1X,A24,I10)') 'eh+ singlet        =',nEhpS
  write(*,'(1X,A24,I10)') 'eh+ triplet        =',nEhpT
  write(*,'(1X,A24,I10)') 'pp+ singlet        =',nEepS
  write(*,'(1X,A24,I10)') 'pp+ triplet        =',nEepT
  write(*,'(1X,A24,I10)') 'Total dimension    =',nH
  write(*,*)

! ======================================================================
! Physical block
! ======================================================================

  do p=1,nOrb
    H(p,p) = eHF(p)
  end do

! ======================================================================
! B- : Part I
! ======================================================================

  ib = oBm

  do i=nC+1,nO
    do a=nO+1,nOrb-nR

      ib = ib + 1
      H(ib,ib) = 2d0*eQP(i) - eQP(a)

      do p=nC+1,nOrb-nR
        val = 0.5d0*ERI(i,i,p,a) + DBmI(i,a,p)
        H(p,ib) = val
        H(ib,p) = val
      end do

    end do
  end do

! ======================================================================
! B- : Part II
! ======================================================================

  ib = oBm + nIm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(i) + eQP(j) - eQP(a)

        do p=nC+1,nOrb-nR
          val = 0.5d0*sqrt(1.5d0)*(ERI(i,j,p,a)-ERI(i,j,a,p)) &
              + DBmII(i,j,a,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! B- : Part III
! ======================================================================

  ib = oBm + nIm + nIIm

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ib = ib + 1
        H(ib,ib) = eQP(i) + eQP(j) - eQP(a)

        do p=nC+1,nOrb-nR
          val = 0.5d0*sqrt(0.5d0)*(ERI(i,j,p,a)+ERI(i,j,a,p)) &
              + DBmIII(i,j,a,p)
          H(p,ib) = val
          H(ib,p) = val
        end do

      end do
    end do
  end do

! ======================================================================
! eh- singlet / triplet
! ======================================================================

  ib = oEhmS
  do j=nC+1,nO
    do n=1,nS
      ib = ib + 1
      H(ib,ib) = eQP(j) - eh_sing_Om(n)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(0.5d0)*eh_sing_rho(p,j,n) + DehmS(j,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

  ib = oEhmT
  do j=nC+1,nO
    do n=1,nS
      ib = ib + 1
      H(ib,ib) = eQP(j) - eh_trip_Om(n)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(1.5d0)*eh_trip_rho(p,j,n) + DehmT(j,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

! ======================================================================
! pp- singlet / triplet (hh response modes)
! ======================================================================

  ib = oHhmS
  do a=nO+1,nOrb-nR
    do n=1,nOOs
      ib = ib + 1
      H(ib,ib) = hh_sing_Om(n) - eQP(a)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(0.5d0)*hh_sing_rho(p,a,n) + DhhmS(a,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

  ib = oHhmT
  do a=nO+1,nOrb-nR
    do n=1,nOOt
      ib = ib + 1
      H(ib,ib) = hh_trip_Om(n) - eQP(a)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(1.5d0)*hh_trip_rho(p,a,n) + DhhmT(a,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

! ======================================================================
! B+ : Part I
! ======================================================================

  ib = oBp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ib = ib + 1
      H(ib,ib) = 2d0*eQP(a) - eQP(i)
      do p=nC+1,nOrb-nR
        val = 0.5d0*ERI(a,a,p,i) + DBpI(i,a,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

! ======================================================================
! B+ : Part II
! ======================================================================

  ib = oBp + nIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1
        H(ib,ib) = eQP(a) + eQP(b) - eQP(i)
        do p=nC+1,nOrb-nR
          val = 0.5d0*sqrt(1.5d0)*(ERI(a,b,p,i)-ERI(a,b,i,p)) &
              + DBpII(i,a,b,p)
          H(p,ib) = val
          H(ib,p) = val
        end do
      end do
    end do
  end do

! ======================================================================
! B+ : Part III
! ======================================================================

  ib = oBp + nIp + nIIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1
        H(ib,ib) = eQP(a) + eQP(b) - eQP(i)
        do p=nC+1,nOrb-nR
          val = 0.5d0*sqrt(0.5d0)*(ERI(a,b,p,i)+ERI(a,b,i,p)) &
              + DBpIII(i,a,b,p)
          H(p,ib) = val
          H(ib,p) = val
        end do
      end do
    end do
  end do

! ======================================================================
! eh+ singlet / triplet
! ======================================================================

  ib = oEhpS
  do b=nO+1,nOrb-nR
    do n=1,nS
      ib = ib + 1
      H(ib,ib) = eQP(b) + eh_sing_Om(n)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(0.5d0)*eh_sing_rho(b,p,n) + DehpS(b,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

  ib = oEhpT
  do b=nO+1,nOrb-nR
    do n=1,nS
      ib = ib + 1
      H(ib,ib) = eQP(b) + eh_trip_Om(n)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(1.5d0)*eh_trip_rho(b,p,n) + DehpT(b,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

! ======================================================================
! pp+ singlet / triplet (ee response modes)
! ======================================================================

  ib = oEepS
  do i=nC+1,nO
    do n=1,nVVs
      ib = ib + 1
      H(ib,ib) = ee_sing_Om(n) - eQP(i)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(0.5d0)*ee_sing_rho(p,i,n) + DeepS(i,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

  ib = oEepT
  do i=nC+1,nO
    do n=1,nVVt
      ib = ib + 1
      H(ib,ib) = ee_trip_Om(n) - eQP(i)
      do p=nC+1,nOrb-nR
        val = 0.5d0*sqrt(1.5d0)*ee_trip_rho(p,i,n) + DeepT(i,n,p)
        H(p,ib) = val
        H(ib,p) = val
      end do
    end do
  end do

! ======================================================================
! B- Part I <-> eh-
! ======================================================================

  ib = oBm
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ib = ib + 1

      jb = oEhmS
      do k=nC+1,nO
        do n=1,nS
          jb = jb + 1
          val = 0d0
          if(i == k) val = -sqrt(2d0)*eh_sing_rho(i,a,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do

      jb = oEhmT
      do k=nC+1,nO
        do n=1,nS
          jb = jb + 1
          val = 0d0
          if(i == k) val = sqrt(6d0)*eh_trip_rho(i,a,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do
    end do
  end do

! ======================================================================
! B- Part II <-> eh-
! ======================================================================

  ib = oBm + nIm
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ib = ib + 1

        jb = oEhmS
        do k=nC+1,nO
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(j == k) val = val + sqrt(3d0)*eh_sing_rho(i,a,n)
            if(i == k) val = val - sqrt(3d0)*eh_sing_rho(j,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do

        jb = oEhmT
        do k=nC+1,nO
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(j == k) val = val + eh_trip_rho(i,a,n)
            if(i == k) val = val - eh_trip_rho(j,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! B- Part III <-> eh-
! ======================================================================

  ib = oBm + nIm + nIIm
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ib = ib + 1

        jb = oEhmS
        do k=nC+1,nO
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(j == k) val = val - eh_sing_rho(i,a,n)
            if(i == k) val = val - eh_sing_rho(j,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do

        jb = oEhmT
        do k=nC+1,nO
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(j == k) val = val + sqrt(3d0)*eh_trip_rho(i,a,n)
            if(i == k) val = val + sqrt(3d0)*eh_trip_rho(j,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! B- <-> pp-
! ======================================================================

! Part I -> singlet pp-

  ib = oBm
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ib = ib + 1
      jb = oHhmS
      do b=nO+1,nOrb-nR
        do n=1,nOOs
          jb = jb + 1
          val = 0d0
          if(a == b) val = -sqrt(2d0)*hh_sing_rho(i,i,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do
    end do
  end do

! Part II -> triplet pp-

  ib = oBm + nIm
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ib = ib + 1
        jb = oHhmT
        do b=nO+1,nOrb-nR
          do n=1,nOOt
            jb = jb + 1
            val = 0d0
            if(a == b) val = -2d0*hh_trip_rho(i,j,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! Part III -> singlet pp-

  ib = oBm + nIm + nIIm
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ib = ib + 1
        jb = oHhmS
        do b=nO+1,nOrb-nR
          do n=1,nOOs
            jb = jb + 1
            val = 0d0
            if(a == b) val = -2d0*hh_sing_rho(i,j,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! B+ Part I <-> eh+
! ======================================================================

  ib = oBp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ib = ib + 1

      jb = oEhpS
      do c=nO+1,nOrb-nR
        do n=1,nS
          jb = jb + 1
          val = 0d0
          if(a == c) val = sqrt(2d0)*eh_sing_rho(i,a,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do

      jb = oEhpT
      do c=nO+1,nOrb-nR
        do n=1,nS
          jb = jb + 1
          val = 0d0
          if(a == c) val = -sqrt(6d0)*eh_trip_rho(i,a,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do
    end do
  end do

! ======================================================================
! B+ Part II <-> eh+
! ======================================================================

  ib = oBp + nIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1

        jb = oEhpS
        do c=nO+1,nOrb-nR
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(a == c) val = val + sqrt(3d0)*eh_sing_rho(i,b,n)
            if(b == c) val = val - sqrt(3d0)*eh_sing_rho(i,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do

        jb = oEhpT
        do c=nO+1,nOrb-nR
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(a == c) val = val + eh_trip_rho(i,b,n)
            if(b == c) val = val - eh_trip_rho(i,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! B+ Part III <-> eh+
! ======================================================================

  ib = oBp + nIp + nIIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1

        jb = oEhpS
        do c=nO+1,nOrb-nR
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(a == c) val = val + eh_sing_rho(i,b,n)
            if(b == c) val = val + eh_sing_rho(i,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do

        jb = oEhpT
        do c=nO+1,nOrb-nR
          do n=1,nS
            jb = jb + 1
            val = 0d0
            if(a == c) val = val - sqrt(3d0)*eh_trip_rho(i,b,n)
            if(b == c) val = val - sqrt(3d0)*eh_trip_rho(i,a,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! B+ <-> pp+
! ======================================================================

! Part I -> singlet pp+

  ib = oBp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ib = ib + 1
      jb = oEepS
      do j=nC+1,nO
        do n=1,nVVs
          jb = jb + 1
          val = 0d0
          if(i == j) val = sqrt(2d0)*ee_sing_rho(a,a,n)
          H(ib,jb) = val
          H(jb,ib) = val
        end do
      end do
    end do
  end do

! Part II -> triplet pp+

  ib = oBp + nIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1
        jb = oEepT
        do j=nC+1,nO
          do n=1,nVVt
            jb = jb + 1
            val = 0d0
            if(i == j) val = 2d0*ee_trip_rho(a,b,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! Part III -> singlet pp+

  ib = oBp + nIp + nIIp
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        ib = ib + 1
        jb = oEepS
        do j=nC+1,nO
          do n=1,nVVs
            jb = jb + 1
            val = 0d0
            if(i == j) val = 2d0*ee_sing_rho(a,b,n)
            H(ib,jb) = val
            H(jb,ib) = val
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! Hermiticity / diagonalization
! ======================================================================

  asym = maxval(abs(H-transpose(H)))
  write(*,'(1X,A34,1X,ES12.4)') 'Max |H-H^T| before diagonalization =',asym

  H = 0.5d0*(H + transpose(H))

  call diagonalize_matrix(nH,H,ePA)

! Spectral weights

  do s=1,nH
    do p=nC+1,nOrb-nR
      Zstate(s) = Zstate(s) + H(p,s)**2
    end do
  end do

  sumZ = sum(Zstate)

  write(*,'(1X,A30,F16.10)') 'Sum of active 1p weights      =',sumZ
  write(*,'(1X,A30,I16)')    'Expected active 1p dimension =',nAct
  write(*,*)

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

  write(*,*)'-------------------------------------------------------'
  write(*,*)'       Full restricted PSD ADC-Parquet spectrum         '
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
! SRG-regularized inverse for internal response denominators
! ======================================================================

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
! Bare removal corrections
! ======================================================================

  double precision function DBmI(ii,aa,pp)

    implicit none
    integer,intent(in) :: ii,aa,pp
    integer            :: nn
    double precision   :: dem

    DBmI = 0d0

    do nn=1,nS
      dem = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      DBmI = DBmI - 0.5d0*eh_sing_rho(aa,ii,nn)*eh_sing_rho(ii,pp,nn)*regularized_inverse(dem)

      dem = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      DBmI = DBmI + 1.5d0*eh_trip_rho(aa,ii,nn)*eh_trip_rho(ii,pp,nn)*regularized_inverse(dem)
    end do

    do nn=1,nVVs
      dem = ee_sing_Om(nn)-2d0*eQP(ii)
      DBmI = DBmI - 0.5d0*ee_sing_rho(ii,ii,nn)*ee_sing_rho(pp,aa,nn)*regularized_inverse(dem)
    end do

  end function DBmI


  double precision function DBmII(ii,jj,aa,pp)

    implicit none
    integer,intent(in) :: ii,jj,aa,pp
    integer            :: nn
    double precision   :: d1,d2,x1,x2

    DBmII = 0d0

    do nn=1,nS

      d1 = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      d2 = eQP(aa)-eQP(jj)+eh_sing_Om(nn)
      x1 = eh_sing_rho(aa,ii,nn)*eh_sing_rho(jj,pp,nn)*regularized_inverse(d1)
      x2 = eh_sing_rho(aa,jj,nn)*eh_sing_rho(ii,pp,nn)*regularized_inverse(d2)
      DBmII = DBmII + 0.25d0*sqrt(6d0)*(x1-x2)

      d1 = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      d2 = eQP(aa)-eQP(jj)+eh_trip_Om(nn)
      x1 = eh_trip_rho(aa,ii,nn)*eh_trip_rho(jj,pp,nn)*regularized_inverse(d1)
      x2 = eh_trip_rho(aa,jj,nn)*eh_trip_rho(ii,pp,nn)*regularized_inverse(d2)
      DBmII = DBmII + 0.25d0*sqrt(6d0)*(x1-x2)

    end do

    do nn=1,nVVt
      d1 = ee_trip_Om(nn)-eQP(ii)-eQP(jj)
      DBmII = DBmII - sqrt(1.5d0)*ee_trip_rho(ii,jj,nn)*ee_trip_rho(pp,aa,nn)*regularized_inverse(d1)
    end do

  end function DBmII


  double precision function DBmIII(ii,jj,aa,pp)

    implicit none
    integer,intent(in) :: ii,jj,aa,pp
    integer            :: nn
    double precision   :: d1,d2,x1,x2

    DBmIII = 0d0

    do nn=1,nS

      d1 = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      d2 = eQP(aa)-eQP(jj)+eh_sing_Om(nn)
      x1 = eh_sing_rho(aa,ii,nn)*eh_sing_rho(jj,pp,nn)*regularized_inverse(d1)
      x2 = eh_sing_rho(aa,jj,nn)*eh_sing_rho(ii,pp,nn)*regularized_inverse(d2)
      DBmIII = DBmIII - 0.25d0*sqrt(2d0)*(x1+x2)

      d1 = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      d2 = eQP(aa)-eQP(jj)+eh_trip_Om(nn)
      x1 = eh_trip_rho(aa,ii,nn)*eh_trip_rho(jj,pp,nn)*regularized_inverse(d1)
      x2 = eh_trip_rho(aa,jj,nn)*eh_trip_rho(ii,pp,nn)*regularized_inverse(d2)
      DBmIII = DBmIII + 0.75d0*sqrt(2d0)*(x1+x2)

    end do

    do nn=1,nVVs
      d1 = ee_sing_Om(nn)-eQP(ii)-eQP(jj)
      DBmIII = DBmIII - sqrt(0.5d0)*ee_sing_rho(ii,jj,nn)*ee_sing_rho(pp,aa,nn)*regularized_inverse(d1)
    end do

  end function DBmIII

! ======================================================================
! Bare addition corrections
! ======================================================================

  double precision function DBpI(ii,aa,pp)

    implicit none
    integer,intent(in) :: ii,aa,pp
    integer            :: nn
    double precision   :: dem

    DBpI = 0d0

    do nn=1,nS
      dem = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      DBpI = DBpI - 0.5d0*eh_sing_rho(aa,ii,nn)*eh_sing_rho(pp,aa,nn)*regularized_inverse(dem)

      dem = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      DBpI = DBpI + 1.5d0*eh_trip_rho(aa,ii,nn)*eh_trip_rho(pp,aa,nn)*regularized_inverse(dem)
    end do

    do nn=1,nOOs
      dem = hh_sing_Om(nn)-2d0*eQP(aa)
      DBpI = DBpI + 0.5d0*hh_sing_rho(aa,aa,nn)*hh_sing_rho(pp,ii,nn)*regularized_inverse(dem)
    end do

  end function DBpI


  double precision function DBpII(ii,aa,bb,pp)

    implicit none
    integer,intent(in) :: ii,aa,bb,pp
    integer            :: nn
    double precision   :: d1,d2,x1,x2

    DBpII = 0d0

    do nn=1,nS

      d1 = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      d2 = eQP(bb)-eQP(ii)+eh_sing_Om(nn)
      x1 = eh_sing_rho(aa,ii,nn)*eh_sing_rho(pp,bb,nn)*regularized_inverse(d1)
      x2 = eh_sing_rho(bb,ii,nn)*eh_sing_rho(pp,aa,nn)*regularized_inverse(d2)
      DBpII = DBpII + 0.25d0*sqrt(6d0)*(x1-x2)

      d1 = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      d2 = eQP(bb)-eQP(ii)+eh_trip_Om(nn)
      x1 = eh_trip_rho(aa,ii,nn)*eh_trip_rho(pp,bb,nn)*regularized_inverse(d1)
      x2 = eh_trip_rho(bb,ii,nn)*eh_trip_rho(pp,aa,nn)*regularized_inverse(d2)
      DBpII = DBpII + 0.25d0*sqrt(6d0)*(x1-x2)

    end do

    do nn=1,nOOt
      d1 = hh_trip_Om(nn)-eQP(aa)-eQP(bb)
      DBpII = DBpII + sqrt(1.5d0)*hh_trip_rho(aa,bb,nn)*hh_trip_rho(pp,ii,nn)*regularized_inverse(d1)
    end do

  end function DBpII


  double precision function DBpIII(ii,aa,bb,pp)

    implicit none
    integer,intent(in) :: ii,aa,bb,pp
    integer            :: nn
    double precision   :: d1,d2,x1,x2

    DBpIII = 0d0

    do nn=1,nS

      d1 = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
      d2 = eQP(bb)-eQP(ii)+eh_sing_Om(nn)
      x1 = eh_sing_rho(aa,ii,nn)*eh_sing_rho(pp,bb,nn)*regularized_inverse(d1)
      x2 = eh_sing_rho(bb,ii,nn)*eh_sing_rho(pp,aa,nn)*regularized_inverse(d2)
      DBpIII = DBpIII - 0.25d0*sqrt(2d0)*(x1+x2)

      d1 = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
      d2 = eQP(bb)-eQP(ii)+eh_trip_Om(nn)
      x1 = eh_trip_rho(aa,ii,nn)*eh_trip_rho(pp,bb,nn)*regularized_inverse(d1)
      x2 = eh_trip_rho(bb,ii,nn)*eh_trip_rho(pp,aa,nn)*regularized_inverse(d2)
      DBpIII = DBpIII + 0.75d0*sqrt(2d0)*(x1+x2)

    end do

    do nn=1,nOOs
      d1 = hh_sing_Om(nn)-eQP(aa)-eQP(bb)
      DBpIII = DBpIII + sqrt(0.5d0)*hh_sing_rho(aa,bb,nn)*hh_sing_rho(pp,ii,nn)*regularized_inverse(d1)
    end do

  end function DBpIII

! ======================================================================
! eh dressed-sector corrections
! ======================================================================

  double precision function DehmS(jj,nn,pp)

    implicit none
    integer,intent(in) :: jj,nn,pp
    integer            :: ii,aa
    double precision   :: dem

    DehmS = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR
        dem = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
        DehmS = DehmS + sqrt(2d0)                                  &
          *(0.5d0*ERI(pp,ii,aa,jj)-ERI(pp,ii,jj,aa))               &
          *eh_sing_rho(aa,ii,nn)*regularized_inverse(dem)
      end do
    end do

  end function DehmS


  double precision function DehmT(jj,nn,pp)

    implicit none
    integer,intent(in) :: jj,nn,pp
    integer            :: ii,aa
    double precision   :: dem

    DehmT = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR
        dem = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
        DehmT = DehmT + sqrt(1.5d0)*ERI(pp,ii,aa,jj)                &
          *eh_trip_rho(aa,ii,nn)*regularized_inverse(dem)
      end do
    end do

  end function DehmT


  double precision function DehpS(bb,nn,pp)

    implicit none
    integer,intent(in) :: bb,nn,pp
    integer            :: ii,aa
    double precision   :: dem

    DehpS = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR
        dem = eQP(aa)-eQP(ii)+eh_sing_Om(nn)
        DehpS = DehpS + sqrt(2d0)                                  &
          *(0.5d0*ERI(pp,aa,ii,bb)-ERI(pp,aa,bb,ii))               &
          *eh_sing_rho(aa,ii,nn)*regularized_inverse(dem)
      end do
    end do

  end function DehpS


  double precision function DehpT(bb,nn,pp)

    implicit none
    integer,intent(in) :: bb,nn,pp
    integer            :: ii,aa
    double precision   :: dem

    DehpT = 0d0

    do ii=nC+1,nO
      do aa=nO+1,nOrb-nR
        dem = eQP(aa)-eQP(ii)+eh_trip_Om(nn)
        DehpT = DehpT + sqrt(1.5d0)*ERI(pp,aa,ii,bb)                &
          *eh_trip_rho(aa,ii,nn)*regularized_inverse(dem)
      end do
    end do

  end function DehpT

! ======================================================================
! pp dressed-sector corrections
! ======================================================================

  double precision function DhhmS(cc,nn,pp)

    implicit none
    integer,intent(in) :: cc,nn,pp
    integer            :: aa,bb
    double precision   :: dem

    DhhmS = 0d0

    do aa=nO+1,nOrb-nR
      do bb=nO+1,nOrb-nR
        dem = hh_sing_Om(nn)-eQP(aa)-eQP(bb)
        DhhmS = DhhmS + sqrt(0.5d0)*ERI(pp,cc,aa,bb)                &
          *hh_sing_rho(aa,bb,nn)*regularized_inverse(dem)
      end do
    end do

  end function DhhmS


  double precision function DhhmT(cc,nn,pp)

    implicit none
    integer,intent(in) :: cc,nn,pp
    integer            :: aa,bb
    double precision   :: dem

    DhhmT = 0d0

    do aa=nO+1,nOrb-nR
      do bb=nO+1,nOrb-nR
        dem = hh_trip_Om(nn)-eQP(aa)-eQP(bb)
        DhhmT = DhhmT + sqrt(1.5d0)*ERI(pp,cc,aa,bb)                &
          *hh_trip_rho(aa,bb,nn)*regularized_inverse(dem)
      end do
    end do

  end function DhhmT


  double precision function DeepS(kk,nn,pp)

    implicit none
    integer,intent(in) :: kk,nn,pp
    integer            :: ii,jj
    double precision   :: dem

    DeepS = 0d0

    do ii=nC+1,nO
      do jj=nC+1,nO
        dem = ee_sing_Om(nn)-eQP(ii)-eQP(jj)
        DeepS = DeepS - sqrt(0.5d0)*ERI(pp,kk,ii,jj)                &
          *ee_sing_rho(ii,jj,nn)*regularized_inverse(dem)
      end do
    end do

  end function DeepS


  double precision function DeepT(kk,nn,pp)

    implicit none
    integer,intent(in) :: kk,nn,pp
    integer            :: ii,jj
    double precision   :: dem

    DeepT = 0d0

    do ii=nC+1,nO
      do jj=nC+1,nO
        dem = ee_trip_Om(nn)-eQP(ii)-eQP(jj)
        DeepT = DeepT - sqrt(1.5d0)*ERI(pp,kk,ii,jj)                &
          *ee_trip_rho(ii,jj,nn)*regularized_inverse(dem)
      end do
    end do

  end function DeepT

end subroutine R_ADC_Parquet
