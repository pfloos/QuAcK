subroutine R_CCD(dotest,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! Spin-adapted closed-shell CCD module

  implicit none

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nOcc,nVir
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: ECC,EcCC
  double precision              :: rcond

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: Loo(:,:)
  double precision,allocatable  :: Lvv(:,:)
  double precision,allocatable  :: Woooo(:,:,:,:)
  double precision,allocatable  :: Wvoov(:,:,:,:)
  double precision,allocatable  :: Wvovo(:,:,:,:)

  double precision,allocatable  :: r(:,:,:,:)
  double precision,allocatable  :: t(:,:,:,:)

  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|   Restricted CCD calculation       |'
  write(*,*)'**************************************'
  write(*,*)

! Active occupied and virtual spaces

  nOcc = nO - nC
  nVir = nV - nR

! Form energy denominator

  allocate(eO(nOcc),eV(nVir))
  allocate(delta_OOVV(nOcc,nOcc,nVir,nVir))

  eO(:) = eHF(nC+1:nO)
  eV(:) = eHF(nO+1:nBas-nR)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

! MP2 guess amplitudes.  t(i,j,a,b) is the alpha-beta (spin-free)
! amplitude; the remaining spin blocks follow from spin symmetry.

  allocate(t(nOcc,nOcc,nVir,nVir))
  call RCCD_guess(nBas,nC,nO,nV,nR,ERI,delta_OOVV,t)

  call RCCD_correlation_energy(nBas,nC,nO,nV,nR,ERI,t,EcMP2)
  EcCC = EcMP2
  ECC  = ERHF + EcCC

! Intermediates and residual

  allocate(Loo(nOcc,nOcc),Lvv(nVir,nVir))
  allocate(Woooo(nOcc,nOcc,nOcc,nOcc))
  allocate(Wvoov(nVir,nOcc,nOcc,nVir))
  allocate(Wvovo(nVir,nOcc,nVir,nOcc))
  allocate(r(nOcc,nOcc,nVir,nVir))

! Memory allocation for DIIS

  allocate(error_diis(nOcc*nOcc*nVir*nVir,max_diis))
  allocate(t_diis(nOcc*nOcc*nVir*nVir,max_diis))

! Initialization

  Conv             = 1d0
  nSCF             = 0
  n_diis           = 0
  t_diis(:,:)      = 0d0
  error_diis(:,:)  = 0d0

!------------------------------------------------------------------------
! Main CCD loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RCCD calculation                                 |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(RCCD)','|','Ec(RCCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    nSCF = nSCF + 1

!   The restricted residual is the alpha-beta block of the spin-orbital
!   CCD residual after analytical spin integration.  In particular, all
!   exchange contractions contain the spin-summed combination 2t-t(exch).

    call form_RCCD_residual(nBas,nC,nO,nV,nR,ERI,delta_OOVV,t, &
                            Loo,Lvv,Woooo,Wvoov,Wvovo,r)

    Conv = maxval(abs(r(:,:,:,:)))

!   Jacobi update

    t(:,:,:,:) = t(:,:,:,:) - r(:,:,:,:)/delta_OOVV(:,:,:,:)

!   DIIS extrapolation of the preconditioned residual

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nOcc*nOcc*nVir*nVir, &
                            nOcc*nOcc*nVir*nVir,n_diis, &
                            error_diis,t_diis,-r/delta_OOVV,t)

    if(abs(rcond) < 1d-15) n_diis = 0

!   Compute correlation energy

    call RCCD_correlation_energy(nBas,nC,nO,nV,nR,ERI,t,EcCC)
    ECC = ERHF + EcCC

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECC+ENuc,'|',EcCC,'|',Conv,'|'

  end do

  write(*,*)'----------------------------------------------------'

! Did it actually converge?

  if(Conv > thresh) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'              RCCD energy                           '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A30,1X,F15.10)')' E(RCCD) = ',ECC+ENuc
  write(*,'(1X,A30,1X,F15.10)')' Ec(RCCD) = ',EcCC
  write(*,*)'----------------------------------------------------'
  write(*,*)

  write(*,'(1X,A15,1X,F10.6)') 'Ec(RMP2) = ',EcMP2
  write(*,*)

! Testing zone.  Keep the historical test label used by QuAcK.

  if(dotest) then

    call dump_test_value('R','CCD correlation energy',EcCC)

  end if

end subroutine R_CCD


subroutine RCCD_guess(nBas,nC,nO,nV,nR,ERI,delta,t)

! Build the restricted MP2 amplitudes.

  implicit none

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: delta(nO-nC,nO-nC,nV-nR,nV-nR)

  integer                       :: i,j,a,b

  double precision,intent(out)  :: t(nO-nC,nO-nC,nV-nR,nV-nR)

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR

          t(i,j,a,b) = -ERI(nC+i,nC+j,nO+a,nO+b)/delta(i,j,a,b)

        end do
      end do
    end do
  end do

end subroutine RCCD_guess


subroutine RCCD_correlation_energy(nBas,nC,nO,nV,nR,ERI,t,EcRCCD)

! Compute the spin-adapted closed-shell CCD correlation energy.

  implicit none

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: t(nO-nC,nO-nC,nV-nR,nV-nR)

  integer                       :: i,j,a,b

  double precision,intent(out)  :: EcRCCD

  EcRCCD = 0d0

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR

          EcRCCD = EcRCCD + ERI(nC+i,nC+j,nO+a,nO+b) &
                            *(2d0*t(i,j,a,b) - t(i,j,b,a))

        end do
      end do
    end do
  end do

end subroutine RCCD_correlation_energy


subroutine form_RCCD_residual(nBas,nC,nO,nV,nR,ERI,delta,t, &
                              Loo,Lvv,Woooo,Wvoov,Wvovo,r)

! Form the spin-adapted closed-shell CCD residual.  These equations are
! equivalent to the alpha-beta block of the antisymmetrized spin-orbital
! equations, but all tensors and contractions use spatial orbitals only.

  implicit none

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: delta(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: t(nO-nC,nO-nC,nV-nR,nV-nR)

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  double precision,intent(out)  :: Loo(nO-nC,nO-nC)
  double precision,intent(out)  :: Lvv(nV-nR,nV-nR)
  double precision,intent(out)  :: Woooo(nO-nC,nO-nC,nO-nC,nO-nC)
  double precision,intent(out)  :: Wvoov(nV-nR,nO-nC,nO-nC,nV-nR)
  double precision,intent(out)  :: Wvovo(nV-nR,nO-nC,nV-nR,nO-nC)
  double precision,intent(out)  :: r(nO-nC,nO-nC,nV-nR,nV-nR)

! Occupied-occupied one-body intermediate

  Loo(:,:) = 0d0

  do k=1,nO-nC
    do i=1,nO-nC
      do l=1,nO-nC
        do c=1,nV-nR
          do d=1,nV-nR

            Loo(k,i) = Loo(k,i) &
              + (2d0*ERI(nC+k,nC+l,nO+c,nO+d) &
                    - ERI(nC+k,nC+l,nO+d,nO+c))*t(i,l,c,d)

          end do
        end do
      end do
    end do
  end do

! Virtual-virtual one-body intermediate

  Lvv(:,:) = 0d0

  do a=1,nV-nR
    do c=1,nV-nR
      do k=1,nO-nC
        do l=1,nO-nC
          do d=1,nV-nR

            Lvv(a,c) = Lvv(a,c) &
              - (2d0*ERI(nC+k,nC+l,nO+c,nO+d) &
                    - ERI(nC+k,nC+l,nO+d,nO+c))*t(k,l,a,d)

          end do
        end do
      end do
    end do
  end do

! Four-occupied intermediate

  do k=1,nO-nC
    do l=1,nO-nC
      do i=1,nO-nC
        do j=1,nO-nC

          Woooo(k,l,i,j) = ERI(nC+k,nC+l,nC+i,nC+j)

          do c=1,nV-nR
            do d=1,nV-nR

              Woooo(k,l,i,j) = Woooo(k,l,i,j) &
                + ERI(nC+k,nC+l,nO+c,nO+d)*t(i,j,c,d)

            end do
          end do

        end do
      end do
    end do
  end do

! Mixed two-body intermediates

  do a=1,nV-nR
    do k=1,nO-nC
      do i=1,nO-nC
        do c=1,nV-nR

          Wvoov(a,k,i,c) = ERI(nC+k,nO+a,nO+c,nC+i)

          do l=1,nO-nC
            do d=1,nV-nR

              Wvoov(a,k,i,c) = Wvoov(a,k,i,c) &
                - 0.5d0*ERI(nC+l,nC+k,nO+d,nO+c)*t(i,l,d,a) &
                - 0.5d0*ERI(nC+l,nC+k,nO+c,nO+d)*t(i,l,a,d) &
                +        ERI(nC+l,nC+k,nO+d,nO+c)*t(i,l,a,d)

            end do
          end do

        end do
      end do
    end do
  end do

  do a=1,nV-nR
    do k=1,nO-nC
      do c=1,nV-nR
        do i=1,nO-nC

          Wvovo(a,k,c,i) = ERI(nC+k,nO+a,nC+i,nO+c)

          do l=1,nO-nC
            do d=1,nV-nR

              Wvovo(a,k,c,i) = Wvovo(a,k,c,i) &
                - 0.5d0*ERI(nC+l,nC+k,nO+c,nO+d)*t(i,l,d,a)

            end do
          end do

        end do
      end do
    end do
  end do

! Residual

  do i=1,nO-nC
    do j=1,nO-nC
      do a=1,nV-nR
        do b=1,nV-nR

          r(i,j,a,b) = ERI(nC+i,nC+j,nO+a,nO+b) &
                     + delta(i,j,a,b)*t(i,j,a,b)

!         Particle-particle ladder

          do c=1,nV-nR
            do d=1,nV-nR

              r(i,j,a,b) = r(i,j,a,b) &
                + ERI(nO+a,nO+b,nO+c,nO+d)*t(i,j,c,d)

            end do
          end do

!         Hole-hole ladder and its quadratic dressing

          do k=1,nO-nC
            do l=1,nO-nC

              r(i,j,a,b) = r(i,j,a,b) + Woooo(k,l,i,j)*t(k,l,a,b)

            end do
          end do

!         One-body insertions

          do c=1,nV-nR

            r(i,j,a,b) = r(i,j,a,b) &
              + Lvv(a,c)*t(i,j,c,b) + Lvv(b,c)*t(j,i,c,a)

          end do

          do k=1,nO-nC

            r(i,j,a,b) = r(i,j,a,b) &
              - Loo(k,i)*t(k,j,a,b) - Loo(k,j)*t(k,i,b,a)

          end do

!         Ring, crossed-ring, and quadratic mixed contractions

          do k=1,nO-nC
            do c=1,nV-nR

              r(i,j,a,b) = r(i,j,a,b) &
                + (2d0*Wvoov(a,k,i,c) - Wvovo(a,k,c,i))*t(k,j,c,b) &
                + (2d0*Wvoov(b,k,j,c) - Wvovo(b,k,c,j))*t(k,i,c,a) &
                - Wvoov(a,k,i,c)*t(k,j,b,c) &
                - Wvoov(b,k,j,c)*t(k,i,a,c) &
                - Wvovo(b,k,c,i)*t(k,j,a,c) &
                - Wvovo(a,k,c,j)*t(k,i,b,c)

            end do
          end do

        end do
      end do
    end do
  end do

end subroutine form_RCCD_residual

