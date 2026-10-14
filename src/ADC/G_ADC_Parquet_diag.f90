subroutine G_ADC_Parquet_diag(nOrb,nC,nO,nV,nR,nS,nOO,nVV,         &
                              eHF,eQP,ERI,eh_Phi,pp_Phi,              &
                              eh_Om,XpY,XmY,eh_rho,                   &
                              ee_Om,X1,Y1,ee_rho,                     &
                              hh_Om,X2,Y2,hh_rho,                     &
                              eADC,ZADC)

! Exact biorthogonal ADC representation of the static-kernel parquet self-energy
! in the generalized (spin-orbital) formalism.
!
! The auxiliary space contains only ADC(3)-class configurations:
!
!   removal : 2h1p (bare) + h x eh-mode + p x hh-mode
!   addition: 2p1h (bare) + p x eh-mode + h x ee-mode
!
! The matrix is non-Hermitian because the exact parquet self-energy has
! different left and right half-amplitudes.  No denominators of the form
! (Omega - Delta) are introduced.
!
! This routine is intended as a diagnostic/prototype analogous to
! G_ADC3_G3W2_diag.f90.  It should be called while XpY/XmY and
! X1/Y1/X2/Y2 are still available in G_evParquet.
!
! The one-body block uses eHF while the auxiliary pole energies use eQP,
! consistently with an upfolding of the self-energy evaluated from the
! current parquet quasiparticle energies.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nOO,nVV
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Phi(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: eh_Om(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  double precision,intent(in)   :: eh_rho(nOrb,nOrb,nS)

  double precision,intent(in)   :: ee_Om(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: ee_rho(nOrb,nOrb,nVV)

  double precision,intent(in)   :: hh_Om(nOO)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)
  double precision,intent(in)   :: hh_rho(nOrb,nOrb,nOO)

! Output variables

  double precision,intent(out)  :: eADC(nOrb)
  double precision,intent(out)  :: ZADC(nOrb)

! Local variables

  integer                       :: nOcc,nVir
  integer                       :: n2h1p,n2p1h
  integer                       :: neh_m,nhh_m,neh_p,nee_p
  integer                       :: nH

  integer                       :: p,i,j,a,b,n,s,idx
  integer                       :: col
  integer                       :: ia,jb,ab,cd,ij,kl
  double precision              :: X,Y
  double precision              :: uL,uR
  double precision              :: cutoff

  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: VL(:,:),VR(:,:)
  double precision,allocatable  :: eval(:),weight(:)

  double precision,allocatable  :: eh_rho0(:,:,:)
  double precision,allocatable  :: ee_rho0(:,:,:)
  double precision,allocatable  :: hh_rho0(:,:,:)

! Useful dimensions

  nOcc = nO - nC
  nVir = nOrb - nR - nO

  n2h1p = nOcc*(nOcc-1)*nVir/2
  n2p1h = nOcc*nVir*(nVir-1)/2

  neh_m = nOcc*nS
  nhh_m = nVir*nOO

  neh_p = nVir*nS
  nee_p = nOcc*nVV

  nH = 1 + n2h1p + neh_m + nhh_m + n2p1h + neh_p + nee_p

  cutoff = 1d-2

! Initialize outputs

  eADC(:) = eHF(:)
  ZADC(:) = 1d0

! Build the bare half-amplitudes directly from the BSE eigenvectors.
! This avoids reconstructing X/Y from rho through dangerous
! (Omega - Delta) denominators.

  allocate(eh_rho0(nOrb,nOrb,nS))
  allocate(ee_rho0(nOrb,nOrb,nVV))
  allocate(hh_rho0(nOrb,nOrb,nOO))

  eh_rho0(:,:,:) = 0d0
  ee_rho0(:,:,:) = 0d0
  hh_rho0(:,:,:) = 0d0

! eh bare half-amplitudes

  do p=nC+1,nOrb-nR
    do i=nC+1,nOrb-nR
      do ia=1,nS

        jb = 0
        do j=nC+1,nO
          do b=nO+1,nOrb-nR
            jb = jb + 1

            X = 0.5d0*(XpY(ia,jb) + XmY(ia,jb))
            Y = 0.5d0*(XpY(ia,jb) - XmY(ia,jb))

            eh_rho0(p,i,ia) = eh_rho0(p,i,ia)                          &
                            + (ERI(p,b,i,j) - ERI(p,b,j,i))*X           &
                            + (ERI(p,j,i,b) - ERI(p,j,b,i))*Y

          end do
        end do

      end do
    end do
  end do

! pp bare half-amplitudes: ee modes

  do p=nC+1,nOrb-nR
    do i=nC+1,nOrb-nR

      do ab=1,nVV

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a+1,nOrb-nR
            cd = cd + 1
            ee_rho0(p,i,ab) = ee_rho0(p,i,ab)                           &
                            + (ERI(p,i,a,b) - ERI(p,i,b,a))*X1(cd,ab)
          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j+1,nO
            kl = kl + 1
            ee_rho0(p,i,ab) = ee_rho0(p,i,ab)                           &
                            + (ERI(p,i,j,b) - ERI(p,i,b,j))*Y1(kl,ab)
          end do
        end do

      end do

! pp bare half-amplitudes: hh modes

      do ij=1,nOO

        cd = 0
        do a=nO+1,nOrb-nR
          do b=a+1,nOrb-nR
            cd = cd + 1
            hh_rho0(p,i,ij) = hh_rho0(p,i,ij)                           &
                            + (ERI(p,i,a,b) - ERI(p,i,b,a))*X2(cd,ij)
          end do
        end do

        kl = 0
        do j=nC+1,nO
          do b=j+1,nO
            kl = kl + 1
            hh_rho0(p,i,ij) = hh_rho0(p,i,ij)                           &
                            + (ERI(p,i,j,b) - ERI(p,i,b,j))*Y2(kl,ij)
          end do
        end do

      end do

    end do
  end do

! Allocate the per-orbital biorthogonal supermatrix

  allocate(H(nH,nH),VL(nH,nH),VR(nH,nH),eval(nH),weight(nH))

  write(*,*)
  write(*,*)'*******************************************'
  write(*,*)'* Generalized exact ADC-Parquet prototype *'
  write(*,*)'*******************************************'
  write(*,*)
  write(*,'(1X,A20,I8)') 'Auxiliary dimension =',nH-1
  write(*,'(1X,A20,I8)') 'Total dimension     =',nH
  write(*,*)

! Main loop over physical orbitals

  do p=nC+1,nOrb-nR

    H(:,:)  = 0d0
    VL(:,:) = 0d0
    VR(:,:) = 0d0
    eval(:) = 0d0

!   ----------------
!   One-body block F
!   ----------------

    H(1,1) = eHF(p)

    col = 1

!   ------------------
!   Bare 2h1p sector
!   ------------------

    do i=nC+1,nO
      do j=i+1,nO
        do a=nO+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(i) + eQP(j) - eQP(a)

!         Left coupling: <pa||ij>

          uL = ERI(p,a,i,j) - ERI(p,a,j,i)

!         Right coupling:
!           <ij||pa>
!         + Gamma_eh(ij,ap) - Gamma_eh(ji,ap)
!         + Gamma_pp(ij,ap)

          uR = (ERI(i,j,p,a) - ERI(i,j,a,p))                            &
             + Gamma_eh(i,j,a,p) - Gamma_eh(j,i,a,p)                   &
             + Gamma_pp(i,j,a,p)

          H(1,col) = uL
          H(col,1) = uR

        end do
      end do
    end do

!   -------------------------
!   eh-dressed removal sector
!   pole = e_j - Omega_eh
!   -------------------------

    do j=nC+1,nO
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(j) - eh_Om(n)

        H(1,col) = eh_rho0(p,j,n)
        H(col,1) = eh_rho (p,j,n)

      end do
    end do

!   -------------------------
!   hh-dressed removal sector
!   pole = Omega_hh - e_a
!   -------------------------

    do a=nO+1,nOrb-nR
      do n=1,nOO

        col = col + 1

        H(col,col) = hh_Om(n) - eQP(a)

        H(1,col) = hh_rho0(p,a,n)
        H(col,1) = hh_rho (p,a,n)

      end do
    end do

!   ------------------
!   Bare 2p1h sector
!   ------------------

    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=a+1,nOrb-nR

          col = col + 1

          H(col,col) = eQP(a) + eQP(b) - eQP(i)

!         Left coupling: <pi||ab>

          uL = ERI(p,i,a,b) - ERI(p,i,b,a)

!         Right coupling:
!           <ab||pi>
!         + Gamma_eh(ab,ip) - Gamma_eh(ba,ip)
!         + Gamma_pp(ab,ip)

          uR = (ERI(a,b,p,i) - ERI(a,b,i,p))                            &
             + Gamma_eh(a,b,i,p) - Gamma_eh(b,a,i,p)                   &
             + Gamma_pp(a,b,i,p)

          H(1,col) = uL
          H(col,1) = uR

        end do
      end do
    end do

!   -------------------------
!   eh-dressed addition sector
!   pole = e_b + Omega_eh
!   -------------------------

    do b=nO+1,nOrb-nR
      do n=1,nS

        col = col + 1

        H(col,col) = eQP(b) + eh_Om(n)

        H(1,col) = eh_rho0(b,p,n)
        H(col,1) = eh_rho (b,p,n)

      end do
    end do

!   -------------------------
!   ee-dressed addition sector
!   pole = Omega_ee - e_i
!   -------------------------

    do i=nC+1,nO
      do n=1,nVV

        col = col + 1

        H(col,col) = ee_Om(n) - eQP(i)

        H(1,col) = ee_rho0(p,i,n)
        H(col,1) = ee_rho (p,i,n)

      end do
    end do

    if(col /= nH) then
      write(*,*) 'G_ADC_Parquet_diag: inconsistent matrix dimension'
      write(*,*) 'col = ',col,' nH = ',nH
      stop
    end if

!   Non-Hermitian diagonalization with left and right eigenvectors.
!   diagonalize_general_matrix_LR bi-normalizes VL^T VR = I
!   for non-degenerate real eigenvalues.

    call diagonalize_general_matrix_LR(nH,H,eval,VL,VR)

!   Signed one-particle spectral weights:
!       Z_s = <L_s|p><p|R_s>
!
!   Negative values are allowed for the original non-PSD parquet
!   self-energy and are therefore a useful diagnostic.

    do s=1,nH
      weight(s) = VL(1,s)*VR(1,s)
    end do

!   Identify the dominant quasiparticle pole by its absolute
!   one-particle weight.

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
  deallocate(eh_rho0,ee_rho0,hh_rho0)

contains

! ======================================================================
! Static generalized parquet kernels
! ======================================================================

  double precision function Gamma_eh(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

!   i Gamma^eh_{p1 p2 p3 p4}
!     = <p1 p2||p3 p4>
!       - Phi^eh_{p1 p2 p4 p3}
!       + Phi^pp_{p1 p2 p3 p4}

    Gamma_eh = ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)                    &
             - eh_Phi(p1,p2,p4,p3)                                    &
             + pp_Phi(p1,p2,p3,p4)

  end function Gamma_eh


  double precision function Gamma_pp(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

!   i Gamma^pp_{p1 p2 p3 p4}
!     = <p1 p2||p3 p4>
!       + Phi^eh_{p1 p2 p3 p4}
!       - Phi^eh_{p1 p2 p4 p3}

    Gamma_pp = ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)                    &
             + eh_Phi(p1,p2,p3,p4)                                    &
             - eh_Phi(p1,p2,p4,p3)

  end function Gamma_pp

end subroutine G_ADC_Parquet_diag

