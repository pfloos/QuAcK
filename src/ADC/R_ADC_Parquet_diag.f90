subroutine R_ADC_Parquet_diag(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                              eHF,eQP,ERI,                                     &
                              eh_sing_Phi,eh_trip_Phi,pp_sing_Phi,pp_trip_Phi, &
                              eh_sing_Om,sing_XpY,sing_XmY,eh_sing_rho,         &
                              eh_trip_Om,trip_XpY,trip_XmY,eh_trip_rho,         &
                              ee_sing_Om,X1s,Y1s,ee_sing_rho,                   &
                              hh_sing_Om,X2s,Y2s,hh_sing_rho,                   &
                              ee_trip_Om,X1t,Y1t,ee_trip_rho,                   &
                              hh_trip_Om,X2t,Y2t,hh_trip_rho,                   &
                              eADC,ZADC)

! Exact spin-adapted, biorthogonal ADC representation of the
! static-kernel restricted parquet self-energy.
!
! This is a diagnostic/prototype analogous to R_ADC3_G3W2, but it
! represents the EXISTING parquet self-energy exactly (in the
! unregularized theory).  No PSD completion is made here.
!
! The matrix is non-Hermitian because the bare and parquet-dressed
! half-amplitudes differ.
!
! Spin weights:
!
!      singlet : wS = 1/2
!      triplet : wT = 3/2
!
! For each eh or pp spin channel we introduce separately
!
!      bare pole     :  - w * g0 * g
!      dressed pole  :  + w * rho0 * rho
!
! This form never divides by Omega - Delta.
!
! The standard GF2 part is represented with the same spatial-orbital
! left/right couplings as R_IPEA_ADC2.f90.
!
! IMPORTANT:
! This routine needs the BSE eigenvectors X/Y.  In R_evParquet it
! should therefore be called after all screened integrals rho have
! been built but BEFORE sing_XpY/sing_XmY, trip_XpY/trip_XmY,
! X1s/Y1s/X2s/Y2s, and X1t/Y1t/X2t/Y2t are deallocated.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  integer,intent(in)            :: nOOs,nVVs,nOOt,nVVt

  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: eh_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_trip_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_sing_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_trip_Phi(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: sing_XpY(nS,nS)
  double precision,intent(in)   :: sing_XmY(nS,nS)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)

  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: trip_XpY(nS,nS)
  double precision,intent(in)   :: trip_XmY(nS,nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)

  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: X1s(nVVs,nVVs)
  double precision,intent(in)   :: Y1s(nOOs,nVVs)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)

  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: X2s(nVVs,nOOs)
  double precision,intent(in)   :: Y2s(nOOs,nOOs)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)

  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: X1t(nVVt,nVVt)
  double precision,intent(in)   :: Y1t(nOOt,nVVt)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)

  double precision,intent(in)   :: hh_trip_Om(nOOt)
  double precision,intent(in)   :: X2t(nVVt,nOOt)
  double precision,intent(in)   :: Y2t(nOOt,nOOt)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOOt)

! Output variables

  double precision,intent(out)  :: eADC(nOrb)
  double precision,intent(out)  :: ZADC(nOrb)

! Local variables

  integer                       :: nOcc,nVir
  integer                       :: nGF2m,nGF2p
  integer                       :: nEHbareM,nEHdressM,nEHbareP,nEHdressP
  integer                       :: nPPbareM,nPPdressM,nPPbareP,nPPdressP
  integer                       :: nH

  integer                       :: p,i,j,a,b,n,s,idx,col
  integer                       :: ia,jb,ab,cd,ij,kl
  double precision              :: X,Y
  double precision              :: wS,wT,sqwS,sqwT
  double precision              :: g0,g
  double precision              :: norm
  double precision              :: cutoff

  double precision,allocatable  :: H(:,:),VL(:,:),VR(:,:)
  double precision,allocatable  :: eval(:),weight(:)

  double precision,allocatable  :: eh_sing_rho0(:,:,:)
  double precision,allocatable  :: eh_trip_rho0(:,:,:)
  double precision,allocatable  :: ee_sing_rho0(:,:,:)
  double precision,allocatable  :: hh_sing_rho0(:,:,:)
  double precision,allocatable  :: ee_trip_rho0(:,:,:)
  double precision,allocatable  :: hh_trip_rho0(:,:,:)

! Constants and dimensions

  wS   = 0.5d0
  wT   = 1.5d0
  sqwS = sqrt(wS)
  sqwT = sqrt(wT)

  cutoff = 1d-2

  nOcc = nO - nC
  nVir = nOrb - nR - nO

! Standard spatial-orbital GF2 sectors

  nGF2m = nOcc*nOcc*nVir
  nGF2p = nOcc*nVir*nVir

! eh channel: for each spin manifold, bare and dressed pole sectors

  nEHbareM  = 2*nOcc*nS
  nEHdressM = 2*nOcc*nS

  nEHbareP  = 2*nVir*nS
  nEHdressP = 2*nVir*nS

! pp channel: singlet + triplet, bare and dressed

  nPPbareM  = nVir*(nOOs+nOOt)
  nPPdressM = nVir*(nOOs+nOOt)

  nPPbareP  = nOcc*(nVVs+nVVt)
  nPPdressP = nOcc*(nVVs+nVVt)

  nH = 1 + nGF2m + nGF2p                                      &
         + nEHbareM + nEHdressM + nEHbareP + nEHdressP        &
         + nPPbareM + nPPdressM + nPPbareP + nPPdressP

! Initialize

  eADC(:) = eHF(:)
  ZADC(:) = 1d0

! ----------------------------------------------------------------------
! Bare transition densities rho0
! ----------------------------------------------------------------------

  allocate(eh_sing_rho0(nOrb,nOrb,nS))
  allocate(eh_trip_rho0(nOrb,nOrb,nS))

  allocate(ee_sing_rho0(nOrb,nOrb,nVVs))
  allocate(hh_sing_rho0(nOrb,nOrb,nOOs))

  allocate(ee_trip_rho0(nOrb,nOrb,nVVt))
  allocate(hh_trip_rho0(nOrb,nOrb,nOOt))

  eh_sing_rho0(:,:,:) = 0d0
  eh_trip_rho0(:,:,:) = 0d0

  ee_sing_rho0(:,:,:) = 0d0
  hh_sing_rho0(:,:,:) = 0d0

  ee_trip_rho0(:,:,:) = 0d0
  hh_trip_rho0(:,:,:) = 0d0

! eh singlet and triplet bare half-amplitudes

  do p=nC+1,nOrb-nR
    do i=nC+1,nOrb-nR

      do ia=1,nS

        jb = 0
        do j=nC+1,nO
          do b=nO+1,nOrb-nR
            jb = jb + 1

            X = 0.5d0*(sing_XpY(ia,jb) + sing_XmY(ia,jb))
            Y = 0.5d0*(sing_XpY(ia,jb) - sing_XmY(ia,jb))

            eh_sing_rho0(p,i,ia) = eh_sing_rho0(p,i,ia)              &
              + (2d0*ERI(p,b,i,j) - ERI(p,b,j,i))*X                 &
              + (2d0*ERI(p,j,i,b) - ERI(p,j,b,i))*Y

            X = 0.5d0*(trip_XpY(ia,jb) + trip_XmY(ia,jb))
            Y = 0.5d0*(trip_XpY(ia,jb) - trip_XmY(ia,jb))

            eh_trip_rho0(p,i,ia) = eh_trip_rho0(p,i,ia)              &
              - ERI(p,b,j,i)*X                                      &
              - ERI(p,j,b,i)*Y

          end do
        end do

      end do

    end do
  end do

! pp singlet bare half-amplitudes

  do p=nC+1,nOrb-nR
    do i=nC+1,nOrb-nR

!     ee singlet modes

      do ab=1,nVVs

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a,nOrb-nR
            cd = cd + 1

            norm = 1d0
            if(a == b) norm = 1d0/sqrt(2d0)

            ee_sing_rho0(p,i,ab) = ee_sing_rho0(p,i,ab)              &
              + norm*(ERI(p,i,a,b) + ERI(p,i,b,a))*X1s(cd,ab)

          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j,nO
            kl = kl + 1

            norm = 1d0
            if(j == b) norm = 1d0/sqrt(2d0)

            ee_sing_rho0(p,i,ab) = ee_sing_rho0(p,i,ab)              &
              + norm*(ERI(p,i,j,b) + ERI(p,i,b,j))*Y1s(kl,ab)

          end do
        end do

      end do

!     hh singlet modes

      do ij=1,nOOs

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a,nOrb-nR
            cd = cd + 1

            norm = 1d0
            if(a == b) norm = 1d0/sqrt(2d0)

            hh_sing_rho0(p,i,ij) = hh_sing_rho0(p,i,ij)              &
              + norm*(ERI(p,i,a,b) + ERI(p,i,b,a))*X2s(cd,ij)

          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j,nO
            kl = kl + 1

            norm = 1d0
            if(j == b) norm = 1d0/sqrt(2d0)

            hh_sing_rho0(p,i,ij) = hh_sing_rho0(p,i,ij)              &
              + norm*(ERI(p,i,j,b) + ERI(p,i,b,j))*Y2s(kl,ij)

          end do
        end do

      end do

    end do
  end do

! pp triplet bare half-amplitudes

  do p=nC+1,nOrb-nR
    do i=nC+1,nOrb-nR

!     ee triplet modes

      do ab=1,nVVt

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a+1,nOrb-nR
            cd = cd + 1

            ee_trip_rho0(p,i,ab) = ee_trip_rho0(p,i,ab)              &
              + (ERI(p,i,a,b) - ERI(p,i,b,a))*X1t(cd,ab)

          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j+1,nO
            kl = kl + 1

            ee_trip_rho0(p,i,ab) = ee_trip_rho0(p,i,ab)              &
              + (ERI(p,i,j,b) - ERI(p,i,b,j))*Y1t(kl,ab)

          end do
        end do

      end do

!     hh triplet modes

      do ij=1,nOOt

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a+1,nOrb-nR
            cd = cd + 1

            hh_trip_rho0(p,i,ij) = hh_trip_rho0(p,i,ij)              &
              + (ERI(p,i,a,b) - ERI(p,i,b,a))*X2t(cd,ij)

          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j+1,nO
            kl = kl + 1

            hh_trip_rho0(p,i,ij) = hh_trip_rho0(p,i,ij)              &
              + (ERI(p,i,j,b) - ERI(p,i,b,j))*Y2t(kl,ij)

          end do
        end do

      end do

    end do
  end do

! ----------------------------------------------------------------------
! Build and diagonalize one non-Hermitian ADC matrix for each p
! ----------------------------------------------------------------------

  allocate(H(nH,nH),VL(nH,nH),VR(nH,nH),eval(nH),weight(nH))

  write(*,*)
  write(*,*)'******************************************'
  write(*,*)'* Restricted exact ADC-Parquet prototype *'
  write(*,*)'******************************************'
  write(*,*)
  write(*,'(1X,A20,I8)') 'Auxiliary dimension =',nH-1
  write(*,'(1X,A20,I8)') 'Total dimension     =',nH
  write(*,*)

  do p=nC+1,nOrb-nR

    H(:,:)      = 0d0
    VL(:,:)     = 0d0
    VR(:,:)     = 0d0
    eval(:)     = 0d0
    weight(:)   = 0d0

    H(1,1) = eHF(p)
    col = 1

!   ================================================================
!   GF2 removal sector: standard R-IPEA-ADC(2) spatial-orbital form
!   ================================================================

    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

          H(1,col) = 2d0*ERI(p,a,i,j) - ERI(p,a,j,i)
          H(col,1) = ERI(p,a,i,j)

        end do
      end do
    end do

!   ================================================================
!   GF2 addition sector
!   ================================================================

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

          H(1,col) = 2d0*ERI(p,i,a,b) - ERI(p,i,b,a)
          H(col,1) = ERI(p,i,a,b)

        end do
      end do
    end do

!   ================================================================
!   eh singlet: bare removal poles
!
!      - 1/2 * g0^S * g^S
!   ================================================================

    do j=nC+1,nO
      do i=nC+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

          g0 = V_eh_sing_X(p,j,i,a)
          g  = G_eh_sing_X(p,j,i,a)

          H(1,col) = +sqwS*g0
          H(col,1) = -sqwS*g

        end do
      end do
    end do

!   ================================================================
!   eh triplet: bare removal poles
!
!      - 3/2 * g0^T * g^T
!   ================================================================

    do j=nC+1,nO
      do i=nC+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

          g0 = V_eh_trip_X(p,j,i,a)
          g  = G_eh_trip_X(p,j,i,a)

          H(1,col) = +sqwT*g0
          H(col,1) = -sqwT*g

        end do
      end do
    end do

!   ================================================================
!   eh singlet: dressed removal poles
!
!      + 1/2 * rho0^S * rho^S
!   ================================================================

    do j=nC+1,nO
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(j) - eh_sing_Om(n)

        H(1,col) = sqwS*eh_sing_rho0(p,j,n)
        H(col,1) = sqwS*eh_sing_rho (p,j,n)

      end do
    end do

!   ================================================================
!   eh triplet: dressed removal poles
!   ================================================================

    do j=nC+1,nO
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(j) - eh_trip_Om(n)

        H(1,col) = sqwT*eh_trip_rho0(p,j,n)
        H(col,1) = sqwT*eh_trip_rho (p,j,n)

      end do
    end do

!   ================================================================
!   eh singlet: bare addition poles
!   ================================================================

    do b=nO+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

          g0 = V_eh_sing_X(b,p,i,a)
          g  = G_eh_sing_X(b,p,i,a)

          H(1,col) = +sqwS*g0
          H(col,1) = -sqwS*g

        end do
      end do
    end do

!   ================================================================
!   eh triplet: bare addition poles
!   ================================================================

    do b=nO+1,nOrb-nR
      do i=nC+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

          g0 = V_eh_trip_X(b,p,i,a)
          g  = G_eh_trip_X(b,p,i,a)

          H(1,col) = +sqwT*g0
          H(col,1) = -sqwT*g

        end do
      end do
    end do

!   ================================================================
!   eh singlet: dressed addition poles
!   ================================================================

    do b=nO+1,nOrb-nR
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(b) + eh_sing_Om(n)

        H(1,col) = sqwS*eh_sing_rho0(b,p,n)
        H(col,1) = sqwS*eh_sing_rho (b,p,n)

      end do
    end do

!   ================================================================
!   eh triplet: dressed addition poles
!   ================================================================

    do b=nO+1,nOrb-nR
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(b) + eh_trip_Om(n)

        H(1,col) = sqwT*eh_trip_rho0(b,p,n)
        H(col,1) = sqwT*eh_trip_rho (b,p,n)

      end do
    end do

!   ================================================================
!   pp singlet: bare hh removal poles
!   ================================================================

    do a=nO+1,nOrb-nR

      ij = 0
      do i=nC+1,nO
        do j=i,nO
          ij = ij + 1

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

          g0 = V_pp_sing(p,a,i,j)
          g  = G_pp_sing(p,a,i,j)

          H(1,col) = +sqwS*g0
          H(col,1) = -sqwS*g

        end do
      end do

    end do

!   ================================================================
!   pp triplet: bare hh removal poles
!   ================================================================

    do a=nO+1,nOrb-nR

      ij = 0
      do i=nC+1,nO
        do j=i+1,nO
          ij = ij + 1

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

          g0 = V_pp_trip(p,a,i,j)
          g  = G_pp_trip(p,a,i,j)

          H(1,col) = +sqwT*g0
          H(col,1) = -sqwT*g

        end do
      end do

    end do

!   ================================================================
!   pp singlet: dressed hh removal poles
!   ================================================================

    do a=nO+1,nOrb-nR
      do n=1,nOOs

        col = col + 1

        H(col,col) = hh_sing_Om(n) - eQP(a)

        H(1,col) = sqwS*hh_sing_rho0(p,a,n)
        H(col,1) = sqwS*hh_sing_rho (p,a,n)

      end do
    end do

!   ================================================================
!   pp triplet: dressed hh removal poles
!   ================================================================

    do a=nO+1,nOrb-nR
      do n=1,nOOt

        col = col + 1

        H(col,col) = hh_trip_Om(n) - eQP(a)

        H(1,col) = sqwT*hh_trip_rho0(p,a,n)
        H(col,1) = sqwT*hh_trip_rho (p,a,n)

      end do
    end do

!   ================================================================
!   pp singlet: bare ee addition poles
!   ================================================================

    do i=nC+1,nO

      ab = 0
      do a=nO+1,nOrb-nR
        do b=a,nOrb-nR
          ab = ab + 1

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

          g0 = V_pp_sing(p,i,a,b)
          g  = G_pp_sing(p,i,a,b)

          H(1,col) = +sqwS*g0
          H(col,1) = -sqwS*g

        end do
      end do

    end do

!   ================================================================
!   pp triplet: bare ee addition poles
!   ================================================================

    do i=nC+1,nO

      ab = 0
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR
          ab = ab + 1

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

          g0 = V_pp_trip(p,i,a,b)
          g  = G_pp_trip(p,i,a,b)

          H(1,col) = +sqwT*g0
          H(col,1) = -sqwT*g

        end do
      end do

    end do

!   ================================================================
!   pp singlet: dressed ee addition poles
!   ================================================================

    do i=nC+1,nO
      do n=1,nVVs

        col = col + 1

        H(col,col) = ee_sing_Om(n) - eQP(i)

        H(1,col) = sqwS*ee_sing_rho0(p,i,n)
        H(col,1) = sqwS*ee_sing_rho (p,i,n)

      end do
    end do

!   ================================================================
!   pp triplet: dressed ee addition poles
!   ================================================================

    do i=nC+1,nO
      do n=1,nVVt

        col = col + 1

        H(col,col) = ee_trip_Om(n) - eQP(i)

        H(1,col) = sqwT*ee_trip_rho0(p,i,n)
        H(col,1) = sqwT*ee_trip_rho (p,i,n)

      end do
    end do

    if(col /= nH) then
      write(*,*) 'R_ADC_Parquet_diag: inconsistent matrix dimension'
      write(*,*) 'col = ',col,' nH = ',nH
      stop
    end if

!   Non-Hermitian diagonalization.
!
!   diagonalize_general_matrix_LR returns left and right eigenvectors
!   with VL^T VR = I for the generic non-degenerate case.

    call diagonalize_general_matrix_LR(nH,H,eval,VL,VR)

!   Signed one-particle spectral weights

    do s=1,nH
      weight(s) = VL(1,s)*VR(1,s)
    end do

!   Dominant quasiparticle pole

    idx = maxloc(abs(weight),dim=1)

    eADC(p) = eval(idx)
    ZADC(p) = weight(idx)

    write(*,*)'---------------------------------------------------------------'
    write(*,'(1X,A18,I4)') 'ADC-Parquet orbital',p
    write(*,'(1X,A15,F16.8,A4)') 'QP energy = ',eADC(p)*HaToeV,' eV'
    write(*,'(1X,A15,F16.8)')    'QP weight = ',ZADC(p)
    write(*,*)' Significant poles:'
    write(*,'(1X,A6,2X,A16,2X,A16)') '#','energy (eV)','signed weight'

    do s=1,nH
      if(abs(weight(s)) > cutoff) then
        write(*,'(1X,I6,2X,F16.8,2X,F16.8)') s,eval(s)*HaToeV,weight(s)
      end if
    end do

    write(*,*)'---------------------------------------------------------------'
    write(*,*)

  end do

  deallocate(H,VL,VR,eval,weight)

  deallocate(eh_sing_rho0,eh_trip_rho0)
  deallocate(ee_sing_rho0,hh_sing_rho0)
  deallocate(ee_trip_rho0,hh_trip_rho0)

contains

! ======================================================================
! Bare and parquet-dressed spin-adapted eh pair couplings
! ======================================================================

  double precision function V_eh_sing_X(p1,q1,i1,a1)

    implicit none
    integer,intent(in) :: p1,q1,i1,a1

    V_eh_sing_X = 2d0*ERI(p1,a1,q1,i1) - ERI(p1,a1,i1,q1)

  end function V_eh_sing_X


  double precision function G_eh_sing_X(p1,q1,i1,a1)

    implicit none
    integer,intent(in) :: p1,q1,i1,a1

    G_eh_sing_X = 2d0*ERI(p1,a1,q1,i1) - ERI(p1,a1,i1,q1)             &
                - 0.5d0*eh_sing_Phi(p1,a1,i1,q1)                     &
                - 1.5d0*eh_trip_Phi(p1,a1,i1,q1)                     &
                + 0.5d0*pp_sing_Phi(p1,a1,q1,i1)                     &
                + 1.5d0*pp_trip_Phi(p1,a1,q1,i1)

  end function G_eh_sing_X


  double precision function V_eh_trip_X(p1,q1,i1,a1)

    implicit none
    integer,intent(in) :: p1,q1,i1,a1

    V_eh_trip_X = -ERI(p1,a1,i1,q1)

  end function V_eh_trip_X


  double precision function G_eh_trip_X(p1,q1,i1,a1)

    implicit none
    integer,intent(in) :: p1,q1,i1,a1

    G_eh_trip_X = -ERI(p1,a1,i1,q1)                                  &
                - 0.5d0*eh_sing_Phi(p1,a1,i1,q1)                     &
                + 0.5d0*eh_trip_Phi(p1,a1,i1,q1)                     &
                - 0.5d0*pp_sing_Phi(p1,a1,q1,i1)                     &
                + 0.5d0*pp_trip_Phi(p1,a1,q1,i1)

  end function G_eh_trip_X

! ======================================================================
! Bare and parquet-dressed spin-adapted pp pair couplings
! ======================================================================

  double precision function V_pp_sing(p1,q1,r1,s1)

    implicit none
    integer,intent(in) :: p1,q1,r1,s1
    double precision   :: fac

    fac = 1d0
    if(r1 == s1) fac = 1d0/sqrt(2d0)

    V_pp_sing = fac*(ERI(p1,q1,r1,s1) + ERI(p1,q1,s1,r1))

  end function V_pp_sing


  double precision function G_pp_sing(p1,q1,r1,s1)

    implicit none
    integer,intent(in) :: p1,q1,r1,s1
    double precision   :: fac

    fac = 1d0
    if(r1 == s1) fac = 1d0/sqrt(2d0)

    G_pp_sing = fac*(                                                     &
          ERI(p1,q1,r1,s1) + ERI(p1,q1,s1,r1)                          &
        + 0.5d0*eh_sing_Phi(p1,q1,r1,s1)                               &
        - 1.5d0*eh_trip_Phi(p1,q1,r1,s1)                               &
        + 0.5d0*eh_sing_Phi(p1,q1,s1,r1)                               &
        - 1.5d0*eh_trip_Phi(p1,q1,s1,r1) )

  end function G_pp_sing


  double precision function V_pp_trip(p1,q1,r1,s1)

    implicit none
    integer,intent(in) :: p1,q1,r1,s1

    V_pp_trip = ERI(p1,q1,r1,s1) - ERI(p1,q1,s1,r1)

  end function V_pp_trip


  double precision function G_pp_trip(p1,q1,r1,s1)

    implicit none
    integer,intent(in) :: p1,q1,r1,s1

    G_pp_trip = ERI(p1,q1,r1,s1) - ERI(p1,q1,s1,r1)                    &
              + 0.5d0*eh_sing_Phi(p1,q1,r1,s1)                         &
              + 0.5d0*eh_trip_Phi(p1,q1,r1,s1)                         &
              - 0.5d0*eh_sing_Phi(p1,q1,s1,r1)                         &
              - 0.5d0*eh_trip_Phi(p1,q1,s1,r1)

  end function G_pp_trip

end subroutine R_ADC_Parquet_diag

