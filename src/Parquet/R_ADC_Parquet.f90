subroutine R_ADC_Parquet(nOrb,nC,nO,nV,nR,eHF,eQP,ERI,                     &
                         eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi, &
                         eADC,ZADC)

! Hermitian / PSD ADC-Parquet in the restricted spin-adapted formalism.
!
! This is a PSD ADC completion of the static parquet kernels, not an
! exact rewriting of the original parquet self-energy.
!
! The auxiliary basis is exactly the three-doublet-sector 2h1p/2p1h
! basis used in QuAcK's full R_IPEA_ADC3 implementation.
!
! The crucial point is that the parquet interaction cannot be represented
! by a single effective spatial Veh/Vpp tensor.  The spin adaptation keeps
! the four independent pair kernels explicitly:
!
!   Gamma_eh_S, Gamma_eh_T, Gamma_pp_S, Gamma_pp_T
!
! At Phi=0 these reduce to the usual singlet/triplet combinations of the
! bare ERIs and every C block reduces exactly to R_IPEA_ADC3.
!
! The full real symmetric matrix is diagonalized once; all one-particle
! spectral weights are therefore non-negative.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR

  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

! Output variables

  double precision,intent(out)  :: eADC(nOrb)
  double precision,intent(out)  :: ZADC(nOrb)

! Local variables

  integer                       :: nOcc,nVir

  integer                       :: nI_2h1p,nII_2h1p,nIII_2h1p
  integer                       :: nI_2p1h,nII_2p1h,nIII_2p1h
  integer                       :: n2h1p,n2p1h,nH,nAct

  integer                       :: p,i,j,k,l,a,b,c,d,s
  integer                       :: ija,klc,iab,kcd
  integer                       :: idx,it

  double precision              :: val
  double precision              :: cutoff
  double precision              :: asym,sumZ
  double precision              :: tol_deg,wgrp,wbest
  double precision              :: start_timing,end_timing,timing

  double precision,external     :: Kronecker_delta

  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: ePA(:)
  double precision,allocatable  :: Zstate(:)

! Dimensions

  nOcc = nO - nC
  nVir = nOrb - nR - nO
  nAct = nOrb - nR - nC

  nI_2h1p   = nOcc*nVir
  nII_2h1p  = nOcc*(nOcc-1)*nVir/2
  nIII_2h1p = nII_2h1p

  nI_2p1h   = nOcc*nVir
  nII_2p1h  = nVir*(nVir-1)*nOcc/2
  nIII_2p1h = nII_2p1h

  n2h1p = nI_2h1p + nII_2h1p + nIII_2h1p
  n2p1h = nI_2p1h + nII_2p1h + nIII_2p1h

  nH = nOrb + n2h1p + n2p1h

  cutoff  = 1d-2
  tol_deg = 1d-8

! Allocation

  allocate(H(nH,nH),ePA(nH),Zstate(nH))

! Initialization

  H(:,:)    = 0d0
  ePA(:)    = 0d0
  Zstate(:) = 0d0

  eADC(:)   = eHF(:)
  ZADC(:)   = 1d0

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Restricted Hermitian ADC-Parquet     *'
  write(*,*)'****************************************'
  write(*,*)
  write(*,'(1X,A20,I8)') '2h1p dimension     =',n2h1p
  write(*,'(1X,A20,I8)') '2p1h dimension     =',n2p1h
  write(*,'(1X,A20,I8)') 'Total dimension    =',nH
  write(*,*)

  call wall_time(start_timing)

! ======================================================================
! F block
! ======================================================================

  do p=1,nOrb
    H(p,p) = eHF(p)
  end do

! ======================================================================
! U^- : Part I
! ======================================================================

  do p=1,nOrb

    ija = nOrb

    do i=nC+1,nO
      do a=nO+1,nOrb-nR

        ija = ija + 1

        val = ERI(i,i,p,a)

        H(p,ija) = val
        H(ija,p) = val

      end do
    end do

  end do

! ======================================================================
! U^- : Part II
! ======================================================================

  do p=1,nOrb

    ija = nOrb + nI_2h1p

    do i=nC+1,nO
      do j=i+1,nO
        do a=nO+1,nOrb-nR

          ija = ija + 1

          val = sqrt(1.5d0)*(ERI(i,j,p,a) - ERI(i,j,a,p))

          H(p,ija) = val
          H(ija,p) = val

        end do
      end do
    end do

  end do

! ======================================================================
! U^- : Part III
! ======================================================================

  do p=1,nOrb

    ija = nOrb + nI_2h1p + nII_2h1p

    do i=nC+1,nO
      do j=i+1,nO
        do a=nO+1,nOrb-nR

          ija = ija + 1

          val = sqrt(0.5d0)*(ERI(i,j,p,a) + ERI(i,j,a,p))

          H(p,ija) = val
          H(ija,p) = val

        end do
      end do
    end do

  end do

! ======================================================================
! U^+ : Part I
! ======================================================================

  do p=1,nOrb

    iab = nOrb + n2h1p

    do i=nC+1,nO
      do a=nO+1,nOrb-nR

        iab = iab + 1

        val = ERI(a,a,p,i)

        H(p,iab) = val
        H(iab,p) = val

      end do
    end do

  end do

! ======================================================================
! U^+ : Part II
! ======================================================================

  do p=1,nOrb

    iab = nOrb + n2h1p + nI_2p1h

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR

          iab = iab + 1

          val = sqrt(1.5d0)*(ERI(a,b,p,i) - ERI(a,b,i,p))

          H(p,iab) = val
          H(iab,p) = val

        end do
      end do
    end do

  end do

! ======================================================================
! U^+ : Part III
! ======================================================================

  do p=1,nOrb

    iab = nOrb + n2h1p + nI_2p1h + nII_2p1h

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR

          iab = iab + 1

          val = sqrt(0.5d0)*(ERI(a,b,p,i) + ERI(a,b,i,p))

          H(p,iab) = val
          H(iab,p) = val

        end do
      end do
    end do

  end do

! ======================================================================
! K^- blocks
! ======================================================================

  ija = nOrb

  do i=nC+1,nO
    do a=nO+1,nOrb-nR

      ija = ija + 1
      H(ija,ija) = 2d0*eQP(i) - eQP(a)

    end do
  end do

  ija = nOrb + nI_2h1p

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ija = ija + 1
        H(ija,ija) = eQP(i) + eQP(j) - eQP(a)

      end do
    end do
  end do

  ija = nOrb + nI_2h1p + nII_2h1p

  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ija = ija + 1
        H(ija,ija) = eQP(i) + eQP(j) - eQP(a)

      end do
    end do
  end do

! ======================================================================
! K^+ blocks
! ======================================================================

  iab = nOrb + n2h1p

  do i=nC+1,nO
    do a=nO+1,nOrb-nR

      iab = iab + 1
      H(iab,iab) = 2d0*eQP(a) - eQP(i)

    end do
  end do

  iab = nOrb + n2h1p + nI_2p1h

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        iab = iab + 1
        H(iab,iab) = eQP(a) + eQP(b) - eQP(i)

      end do
    end do
  end do

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h

  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        iab = iab + 1
        H(iab,iab) = eQP(a) + eQP(b) - eQP(i)

      end do
    end do
  end do

! ======================================================================
! C^- : Part I-I
!
! Bare limit:
!   -delta_ac (ii|kk) - delta_ik[(ci|ka)-2(ci|ak)]
! = -1/2 delta_ac Gamma_pp_S(ii,kk)
!   +delta_ik Gamma_eh_S(ci,ak)
! ======================================================================

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1

      klc = nOrb
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          klc = klc + 1

          H(ija,klc) = H(ija,klc)                                    &
            - 0.5d0*Kronecker_delta(a,c)*Gamma_pp_S(i,i,k,k)         &
            + Kronecker_delta(i,k)*Gamma_eh_S(c,i,a,k)

        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part II-II
! ======================================================================

  ija = nOrb + nI_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        klc = nOrb + nI_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1

              H(ija,klc) = H(ija,klc)                                      &
                - Kronecker_delta(a,c)*Gamma_pp_T(i,j,k,l)                  &
                + Kronecker_delta(j,l)*Eh_mix(c,i,a,k)                      &
                - Kronecker_delta(j,k)*Eh_mix(c,i,a,l)                      &
                + Kronecker_delta(i,k)*Eh_mix(c,j,a,l)                      &
                - Kronecker_delta(i,l)*Eh_mix(c,j,a,k)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part III-III
! ======================================================================

  ija = nOrb + nI_2h1p + nII_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        klc = nOrb + nI_2h1p + nII_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1

              H(ija,klc) = H(ija,klc)                                      &
                - Kronecker_delta(a,c)*Gamma_pp_S(i,j,k,l)                  &
                + 0.5d0*Kronecker_delta(i,k)*Gamma_eh_S(c,j,a,l)            &
                + 0.5d0*Kronecker_delta(j,l)*Gamma_eh_S(c,i,a,k)            &
                + 0.5d0*Kronecker_delta(i,l)*Gamma_eh_S(c,j,a,k)            &
                + 0.5d0*Kronecker_delta(j,k)*Gamma_eh_S(c,i,a,l)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part I-II and II-I
! ======================================================================

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1

      klc = nOrb + nI_2h1p
      do k=nC+1,nO
        do l=k+1,nO
          do c=nO+1,nOrb-nR
            klc = klc + 1

            val = +sqrt(1.5d0)*Kronecker_delta(i,k)*Gamma_eh_T(c,i,a,l) &
                  -sqrt(1.5d0)*Kronecker_delta(i,l)*Gamma_eh_T(c,i,a,k)

            H(ija,klc) = H(ija,klc) + val
            H(klc,ija) = H(klc,ija) + val

          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part I-III and III-I
! ======================================================================

  ija = nOrb
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      ija = ija + 1

      klc = nOrb + nI_2h1p + nII_2h1p
      do k=nC+1,nO
        do l=k+1,nO
          do c=nO+1,nOrb-nR
            klc = klc + 1

            val = ( - Kronecker_delta(a,c)*Gamma_pp_S(i,i,k,l)        &
                    + Kronecker_delta(i,k)*Gamma_eh_S(c,i,a,l)        &
                    + Kronecker_delta(i,l)*Gamma_eh_S(c,i,a,k) )      &
                  /sqrt(2d0)

            H(ija,klc) = H(ija,klc) + val
            H(klc,ija) = H(klc,ija) + val

          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part II-III
! ======================================================================

  ija = nOrb + nI_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        klc = nOrb + nI_2h1p + nII_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1

              H(ija,klc) = H(ija,klc)                                  &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(i,k)*Gamma_eh_T(c,j,a,l) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(j,l)*Gamma_eh_T(c,i,a,k) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(i,l)*Gamma_eh_T(c,j,a,k) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(j,k)*Gamma_eh_T(c,i,a,l)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^- : Part III-II
! ======================================================================

  ija = nOrb + nI_2h1p + nII_2h1p
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        klc = nOrb + nI_2h1p
        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1

              H(ija,klc) = H(ija,klc)                                  &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(i,k)*Gamma_eh_T(c,j,a,l) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(j,l)*Gamma_eh_T(c,i,a,k) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(i,l)*Gamma_eh_T(c,j,a,k) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(j,k)*Gamma_eh_T(c,i,a,l)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part I-I
! ======================================================================

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1

      kcd = nOrb + n2h1p
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          kcd = kcd + 1

          H(iab,kcd) = H(iab,kcd)                                    &
            + 0.5d0*Kronecker_delta(i,k)*Gamma_pp_S(a,a,c,c)         &
            - Kronecker_delta(a,c)*Gamma_eh_S(k,a,i,c)

        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part II-II
! ======================================================================

  iab = nOrb + n2h1p + nI_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1

        kcd = nOrb + n2h1p + nI_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1

              H(iab,kcd) = H(iab,kcd)                                      &
                + Kronecker_delta(i,k)*Gamma_pp_T(a,b,c,d)                  &
                - Kronecker_delta(b,d)*Eh_mix(k,a,i,c)                      &
                + Kronecker_delta(b,c)*Eh_mix(k,a,i,d)                      &
                - Kronecker_delta(a,c)*Eh_mix(k,b,i,d)                      &
                + Kronecker_delta(a,d)*Eh_mix(k,b,i,c)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part III-III
! ======================================================================

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1

        kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1

              H(iab,kcd) = H(iab,kcd)                                      &
                + Kronecker_delta(i,k)*Gamma_pp_S(a,b,c,d)                  &
                - 0.5d0*Kronecker_delta(a,c)*Gamma_eh_S(k,b,i,d)            &
                - 0.5d0*Kronecker_delta(b,d)*Gamma_eh_S(k,a,i,c)            &
                - 0.5d0*Kronecker_delta(a,d)*Gamma_eh_S(k,b,i,c)            &
                - 0.5d0*Kronecker_delta(b,c)*Gamma_eh_S(k,a,i,d)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part I-II and II-I
! ======================================================================

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1

      kcd = nOrb + n2h1p + nI_2p1h
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          do d=c+1,nOrb-nR
            kcd = kcd + 1

            val = -sqrt(1.5d0)*Kronecker_delta(a,c)*Gamma_eh_T(k,a,i,d) &
                  +sqrt(1.5d0)*Kronecker_delta(a,d)*Gamma_eh_T(k,a,i,c)

            H(iab,kcd) = H(iab,kcd) + val
            H(kcd,iab) = H(kcd,iab) + val

          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part I-III and III-I
! ======================================================================

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      iab = iab + 1

      kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
      do k=nC+1,nO
        do c=nO+1,nOrb-nR
          do d=c+1,nOrb-nR
            kcd = kcd + 1

            val = ( + Kronecker_delta(i,k)*Gamma_pp_S(a,a,c,d)        &
                    - Kronecker_delta(a,c)*Gamma_eh_S(k,a,i,d)        &
                    - Kronecker_delta(a,d)*Gamma_eh_S(k,a,i,c) )      &
                  /sqrt(2d0)

            H(iab,kcd) = H(iab,kcd) + val
            H(kcd,iab) = H(kcd,iab) + val

          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part II-III
! ======================================================================

  iab = nOrb + n2h1p + nI_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1

        kcd = nOrb + n2h1p + nI_2p1h + nII_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1

              H(iab,kcd) = H(iab,kcd)                                  &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(a,c)*Gamma_eh_T(k,b,i,d) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(b,d)*Gamma_eh_T(k,a,i,c) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(a,d)*Gamma_eh_T(k,b,i,c) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(b,c)*Gamma_eh_T(k,a,i,d)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! C^+ : Part III-II
! ======================================================================

  iab = nOrb + n2h1p + nI_2p1h + nII_2p1h
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR
        iab = iab + 1

        kcd = nOrb + n2h1p + nI_2p1h
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR
              kcd = kcd + 1

              H(iab,kcd) = H(iab,kcd)                                  &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(a,c)*Gamma_eh_T(k,b,i,d) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(b,d)*Gamma_eh_T(k,a,i,c) &
                + 0.5d0*sqrt(3d0)*Kronecker_delta(a,d)*Gamma_eh_T(k,b,i,c) &
                - 0.5d0*sqrt(3d0)*Kronecker_delta(b,c)*Gamma_eh_T(k,a,i,d)

            end do
          end do
        end do
      end do
    end do
  end do

! ======================================================================
! Enforce exact Hermiticity
! ======================================================================

  asym = maxval(abs(H(:,:) - transpose(H(:,:))))
  write(*,'(1X,A32,1X,ES12.4)') 'Max pre-symmetrization |H-H^T| =',asym

  H(:,:) = 0.5d0*(H(:,:) + transpose(H(:,:)))

  call wall_time(end_timing)
  timing = end_timing - start_timing

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for ADC-PA matrix construction =',timing,' seconds'
  write(*,*)

! ======================================================================
! Diagonalization
! ======================================================================

  call wall_time(start_timing)

  call diagonalize_matrix(nH,H,ePA)

  call wall_time(end_timing)
  timing = end_timing - start_timing

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for ADC-PA diagonalization =',timing,' seconds'
  write(*,*)

! ======================================================================
! Spectral weights
! ======================================================================

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
        if(abs(ePA(it+1) - ePA(s)) > tol_deg) exit
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

! ======================================================================
! Print spectrum
! ======================================================================

  write(*,*)'-------------------------------------------------------'
  write(*,*)'          Restricted PSD ADC-Parquet spectrum           '
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
! Spin-adapted parquet pair kernels.
!
! These are exactly the combinations used by the restricted parquet
! screened-integral routines.  Keeping S and T separately is essential:
! a general parquet vertex cannot be represented by one effective spatial
! two-electron integral.
! ======================================================================

  double precision function Gamma_eh_S(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_eh_S = 2d0*ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)           &
      - 0.5d0*eh_sing_Phi(p1,p2,p4,p3)                             &
      - 1.5d0*eh_trip_Phi(p1,p2,p4,p3)                             &
      + 0.5d0*pp_sing_Phi(p1,p2,p3,p4)                             &
      + 1.5d0*pp_trip_Phi(p1,p2,p3,p4)

  end function Gamma_eh_S


  double precision function Gamma_eh_T(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_eh_T = - ERI(p1,p2,p4,p3)                                &
      - 0.5d0*eh_sing_Phi(p1,p2,p4,p3)                             &
      + 0.5d0*eh_trip_Phi(p1,p2,p4,p3)                             &
      - 0.5d0*pp_sing_Phi(p1,p2,p3,p4)                             &
      + 0.5d0*pp_trip_Phi(p1,p2,p3,p4)

  end function Gamma_eh_T


  double precision function Gamma_pp_S(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_pp_S = ERI(p1,p2,p3,p4) + ERI(p1,p2,p4,p3)              &
      + 0.5d0*eh_sing_Phi(p1,p2,p3,p4)                             &
      - 1.5d0*eh_trip_Phi(p1,p2,p3,p4)                             &
      + 0.5d0*eh_sing_Phi(p1,p2,p4,p3)                             &
      - 1.5d0*eh_trip_Phi(p1,p2,p4,p3)

  end function Gamma_pp_S


  double precision function Gamma_pp_T(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_pp_T = ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)              &
      + 0.5d0*eh_sing_Phi(p1,p2,p3,p4)                             &
      + 0.5d0*eh_trip_Phi(p1,p2,p3,p4)                             &
      - 0.5d0*eh_sing_Phi(p1,p2,p4,p3)                             &
      - 0.5d0*eh_trip_Phi(p1,p2,p4,p3)

  end function Gamma_pp_T


  double precision function Eh_mix(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

!   In the bare limit:
!     1/2 Gamma_eh_S + Gamma_eh_T
!       = (pq|rs) - 3/2 (pq|sr)

    Eh_mix = 0.5d0*Gamma_eh_S(p1,p2,p3,p4) + Gamma_eh_T(p1,p2,p3,p4)

  end function Eh_mix

end subroutine R_ADC_Parquet
