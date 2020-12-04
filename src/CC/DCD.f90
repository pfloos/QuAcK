subroutine DCD(maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! DCD module

  implicit none

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2,EcMP3,EcMP4
  double precision              :: EDCD,EcDCD

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVOV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)

  double precision,allocatable  :: xO(:,:)
  double precision,allocatable  :: xV(:,:)
  double precision,allocatable  :: xOV(:,:,:,:)

  double precision,allocatable  :: r(:,:,:,:)
  double precision,allocatable  :: t(:,:,:,:)
  double precision,allocatable  :: tt(:,:,:,:)

  integer                       :: n_diis
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          DCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Form energy denominator

  allocate(eO(nO-nC),eV(nV-nR))
  allocate(delta(nO-nC,nO-nC,nV-nR,nV-nR))

  eO(:) = eHF(nC+1:nO)
  eV(:) = eHF(nO+1:nBas-nR)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta)

! Create integral batches

  allocate(OOOO(nO-nC,nO-nC,nO-nC,nO-nC), &
           OOVV(nO-nC,nO-nC,nV-nR,nV-nR), & 
           OVOV(nO-nC,nV-nR,nO-nC,nV-nR), & 
           OVVO(nO-nC,nV-nR,nV-nR,nO-nC), &
           VVVV(nV-nR,nV-nR,nV-nR,nV-nR))

  OOOO(:,:,:,:) = ERI(nC+1:nO     ,nC+1:nO     ,nC+1:nO     ,nC+1:nO     )
  OOVV(:,:,:,:) = ERI(nC+1:nO     ,nC+1:nO     ,nO+1:nBas-nR,nO+1:nBas-nR)
  OVOV(:,:,:,:) = ERI(nC+1:nO     ,nO+1:nBas-nR,nC+1:nO     ,nO+1:nBas-nR)
  OVVO(:,:,:,:) = ERI(nC+1:nO     ,nO+1:nBas-nR,nO+1:nBas-nR,nC+1:nO     )
  VVVV(:,:,:,:) = ERI(nO+1:nBas-nR,nO+1:nBas-nR,nO+1:nBas-nR,nO+1:nBas-nR)
 
! MP2 guess amplitudes

  allocate(t(nO-nC,nO-nC,nV-nR,nV-nR))

  t(:,:,:,:) = -OOVV(:,:,:,:)/delta(:,:,:,:)

! Memory allocation for DIIS

  allocate(error_diis((nO-nR)**2*(nV-nR)**2,max_diis), & 
               t_diis((nO-nR)**2*(nV-nR)**2,max_diis))

! Initialization

  allocate(tt(nO-nC,nO-nC,nV-nR,nV-nR),     &
           r(nO-nC,nO-nC,nV-nR,nV-nR),      &
           xO(nO-nC,nO-nC),xV(nV-nR,nV-nR), &
           xOV(nO-nC,nO-nC,nV-nR,nV-nR))

  Conv = 1d0
  nSCF = 0

  n_diis          = 0
  t_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| DCD calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(DCD)','|','Ec(DCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Compute correlation energy

    EcDCD = 0d0
    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
 
            EcDCD = EcDCD + (2d0*OOVV(i,j,a,b) - OOVV(i,j,b,a))*t(i,j,a,b)
 
          enddo
        enddo
      enddo
    enddo

    ! Dump results

    EDCD = ERHF + EcDCD

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',EDCD+ENuc,'|',EcDCD,'|',Conv,'|'

    ! Increment 

    nSCF = nSCF + 1

    ! Form interemediate arrays

    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
            tt(i,j,a,b) = 2d0*t(i,j,a,b) - t(i,j,b,a)
          enddo
        enddo
      enddo
    enddo

    xV(:,:) = 0d0
    do a=1,nV-nR
      do b=1,nV-nR
        do k=1,nO-nC
          do l=1,nO-nC
            do c=1,nV-nR
              xV(a,b) = xV(a,b) + tt(k,l,a,c)*OOVV(l,k,c,b)
            end do
          end do
        end do
      end do
    end do

    xO(:,:) = 0d0
    do i=1,nO-nC
      do j=1,nO-nC
        do k=1,nO-nC
          do c=1,nV-nR
            do d=1,nV-nR
              xO(i,j) = xO(i,j) + tt(i,k,c,d)*OOVV(k,j,d,c)
            end do
          end do
        end do
      end do
    end do

    xOV(:,:,:,:) = 0d0
    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
            do k=1,nO-nC
              do c=1,nV-nR
                xOV(i,j,a,b) = xOV(i,j,a,b) + tt(i,k,a,c)*OOVV(k,j,c,b)
              end do
            end do
          end do
        end do
      end do
    end do

    ! Compute residual

    r(:,:,:,:) = OOVV(:,:,:,:) + delta(:,:,:,:)*t(:,:,:,:) 

    do i=1,nO-nC
      do j=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR

            do c=1,nV-nR
              r(i,j,a,b) = r(i,j,a,b) - 0.5d0*xV(a,c)*t(i,j,c,b) &
                                      - 0.5d0*xV(b,c)*t(j,i,c,a)
              do d=1,nV-nR
                r(i,j,a,b) = r(i,j,a,b) + VVVV(a,b,c,d)*t(i,j,c,d) 
              end do
            end do

            do k=1,nO-nC
              r(i,j,a,b) = r(i,j,a,b) - 0.5d0*xO(k,i)*t(k,j,a,b) &
                                      - 0.5d0*xO(k,j)*t(k,i,b,a)
              do l=1,nO-nC
                r(i,j,a,b) = r(i,j,a,b) + OOOO(i,j,k,l)*t(k,l,a,b) 
              end do
            end do

            do k=1,nO-nC
              do c=1,nV-nR

                r(i,j,a,b) = r(i,j,a,b) + xOV(i,k,a,c)*tt(k,j,c,b)  &
                                        - OVOV(k,a,i,c)*t(k,j,c,b)  &
                                        - OVOV(k,b,i,c)*t(k,j,a,c)  &
                                        - OVOV(k,b,j,c)*t(k,i,c,a)  &
                                        - OVOV(k,a,j,c)*t(k,i,b,c)  &        
                                        + tt(i,k,a,c)*OVVO(k,b,c,j) &
                                        + tt(j,k,b,c)*OVVO(k,a,c,i)
              end do
            end do

          end do
        end do
      end do
    end do

    ! Check convergence 

    Conv = maxval(abs(r(:,:,:,:)))
  
    ! Update amplitudes

    t(:,:,:,:) = t(:,:,:,:) - r(:,:,:,:)/delta(:,:,:,:)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,(nO-nC)**2*(nV-nR)**2,(nO-nC)**2*(nV-nR)**2,n_diis,error_diis,t_diis,-r/delta,t)

    !  Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

  enddo
  write(*,*)'----------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

end subroutine DCD
