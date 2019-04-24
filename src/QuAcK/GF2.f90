subroutine GF2(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,V,e0)

! Perform second-order Green function calculation in diagonal approximation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e0(nBas),V(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nSCF,n_diis
  double precision              :: eps
  double precision              :: Conv
  double precision              :: rcond
  double precision,allocatable  :: eGF2(:),eOld(:),Bpp(:,:,:),error_diis(:,:),e_diis(:,:)

  integer                       :: i,j,a,b,p,q

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Second-order Green function calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(Bpp(nBas,nBas,2),eGF2(nBas),eOld(nBas), &
           error_diis(nBas,max_diis),e_diis(nBas,max_diis))

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

    Bpp(:,:,:) = 0d0

    do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nBas-nR

          eps = eGF2(p) + e0(a) - e0(i) - e0(j)

          Bpp(p,q,1) = Bpp(p,q,1) &
                   + (2d0*V(p,a,i,j) - V(p,a,j,i))*V(q,a,i,j)/eps

        enddo
      enddo
      enddo
    enddo
    enddo

    do p=nC+1,nBas-nR
    do q=nC+1,nBas-nR
      do i=nC+1,nO
        do a=nO+1,nBas-nR
        do b=nO+1,nBas-nR

          eps = eGF2(p) + e0(i) - e0(a) - e0(b)

          Bpp(p,q,2) = Bpp(p,q,2) &
                   + (2d0*V(p,i,a,b) - V(p,i,b,a))*V(q,i,a,b)/eps

        enddo
        enddo
      enddo
    enddo
    enddo

    print*,'Sig2 in GF2'
    call matout(nBas,nBas,Bpp(:,:,1) + Bpp(:,:,2))

!   eGF2(:) = e0(:) &
!           + Bpp(:,1) + Bpp(:,2)

    Conv = maxval(abs(eGF2 - eOld))

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGF2-eOld,eGF2)

!    Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

    eOld = eGF2

    ! Print results

    call print_GF2(nBas,nO,nSCF,Conv,e0,eGF2)

    ! Increment

    nSCF = nSCF + 1

  enddo
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

  endif

end subroutine GF2
