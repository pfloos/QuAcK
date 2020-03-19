subroutine evGF2(maxSCF,thresh,max_diis,linearize,nBas,nC,nO,nV,nR,V,e0)

! Perform eigenvalue self-consistent second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  logical,intent(in)            :: linearize
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: e0(nBas)
  double precision,intent(in)   :: V(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: eps
  double precision              :: Conv
  double precision              :: rcond
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)

  integer                       :: i,j,a,b,p

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Second-order Green function calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(Bpp(nBas,2),Z(nBas),eGF2(nBas),eOld(nBas),error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  Conv            = 1d0
  nSCF            = 0
  n_diis          = 0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGF2(:)         = e0(:)
  eOld(:)         = e0(:)

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Frequency-dependent second-order contribution

    Bpp(:,:) = 0d0

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR

            eps = eGF2(p) + e0(a) - e0(i) - e0(j)

            Bpp(p,1) = Bpp(p,1) &
                     + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(p,a,i,j)/eps

          end do
        end do
      end do
    end do

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR

            eps = eGF2(p) + e0(i) - e0(a) - e0(b)

            Bpp(p,2) = Bpp(p,2) &
                     + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(p,i,a,b)/eps

          end do
        end do
      end do
    end do

    ! Compute the renormalization factor

    Z(:) = 0d0

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nBas-nR

            eps = eGF2(p) + e0(a) - e0(i) - e0(j)

            Z(p) = Z(p) - (2d0*V(p,a,i,j) - V(p,a,j,i))*V(p,a,i,j)/eps**2

          end do
        end do
      end do
    end do

    do p=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          do b=nO+1,nBas-nR

            eps = eGF2(p) + e0(i) - e0(a) - e0(b)

            Z(p) = Z(p) - (2d0*V(p,i,a,b) - V(p,i,b,a))*V(p,i,a,b)/eps**2

          end do
        end do
      end do
    end do

    Z(:) = 1d0/(1d0 - Z(:))

    if(linearize) then

      eGF2(:) = e0(:) + Z(:)*(Bpp(:,1) + Bpp(:,2))

    else

      eGF2(:) = e0(:) + Bpp(:,1) + Bpp(:,2)

    end if

    Conv = maxval(abs(eGF2 - eOld))

    ! Print results

    call print_evGF2(nBas,nO,nSCF,Conv,e0,eGF2)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGF2-eOld,eGF2)

    if(abs(rcond) < 1d-15) n_diis = 0

    eOld(:) = eGF2(:)

    ! Increment

    nSCF = nSCF + 1

  end do
!------------------------------------------------------------------------
! End main SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

  end if

end subroutine evGF2
