subroutine CCGW(maxSCF,thresh,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,e)

! CC-based GW module

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p,q
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  integer                       :: nSCF
  double precision              :: Conv

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_2h1p(:,:,:,:)
  double precision,allocatable  :: delta_2p1h(:,:,:,:)

  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VOOV(:,:,:,:)
  double precision,allocatable  :: NVOO(:,:,:,:)
  double precision,allocatable  :: NOVV(:,:,:,:)

  double precision,allocatable  :: r_2h1p(:,:,:,:)
  double precision,allocatable  :: r_2p1h(:,:,:,:)
  double precision,allocatable  :: t_2h1p(:,:,:,:)
  double precision,allocatable  :: t_2p1h(:,:,:,:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     ring CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Create integral batches

  allocate(OVVO(nO,nV,nV,nO),VOOV(nV,nO,nO,nV),NVOO(nBas,nV,nO,nO),NOVV(nBas,nO,nV,nV))

  OVVO(:,:,:,:) = ERI(   1:nO   ,nO+1:nBas,nO+1:nBas,   1:nO  )
  VOOV(:,:,:,:) = ERI(nO+1:nBas ,   1:nO  ,   1:nO  ,nO+1:nBas)
  NVOO(:,:,:,:) = ERI(   1:nBas ,nO+1:nBas,   1:nO  ,   1:nO  )
  NOVV(:,:,:,:) = ERI(   1:nBas ,   1:nO  ,nO+1:nBas,nO+1:nBas)
 
! Form energy denominator and guess amplitudes

  allocate(eO(nO),eV(nV))
  allocate(delta_2h1p(nO,nO,nV,nBas),delta_2p1h(nO,nV,nV,nBas))
  allocate(t_2h1p(nO,nO,nV,nBas),t_2p1h(nO,nV,nV,nBas))

  eO(:) = e(1:nO)
  eV(:) = e(nO+1:nBas)

  do i=nC+1,nO
    do j=nC+1,nO

      do a=1,nV-nR
        do p=nC+1,nBas-nR

          delta_2h1p(i,j,a,p) = eO(i) + eO(j) - eV(a) - e(p)
          t_2h1p(i,j,a,p) = - sqrt(2d0)*NVOO(p,a,i,j)/delta_2h1p(i,j,a,p)

        end do
      end do

    end do
  end do

  do a=1,nV-nR
    do b=1,nV-nR

      do i=nC+1,nO
        do p=nC+1,nBas-nR

          delta_2p1h(i,a,b,p) = eV(a) + eV(b) - eO(i) - e(p)
          t_2p1h(i,a,b,p) = - sqrt(2d0)*NOVV(p,i,b,a)/delta_2p1h(i,a,b,p)

        end do
      end do

    end do
  end do

! Initialization

  allocate(r_2h1p(nO,nO,nV,nBas),r_2p1h(nO,nV,nV,nBas))
  allocate(eGW(nBas),Z(nBas))

  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| CCGW calculation                                 |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HOMO','|','LUMO','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Compute residual for 2h1p sector

    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR
          do p=nC+1,nBas-nR
 
            r_2h1p(i,j,a,p) = sqrt(2d0)*NVOO(p,a,i,j) + delta_2h1p(i,j,a,p)*t_2h1p(i,j,a,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - 2d0*OVVO(j,c,a,k)*t_2h1p(i,k,c,p)

                do l=nC+1,nO
                  do q=nC+1,nBas-nR

                    r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - sqrt(2d0)*t_2h1p(i,j,a,q)*NVOO(q,c,k,l)*t_2h1p(k,l,c,p)

                  end do
                end do

                do d=1,nV-nR
                  do q=nC+1,nBas-nR

                    r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - sqrt(2d0)*t_2h1p(i,j,a,q)*NOVV(q,k,d,c)*t_2p1h(k,c,d,p)

                  end do
                end do

              end do
            end do
 
          end do
        end do
      end do
    end do
 
!   Compute residual for 2p1h sector
 
    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
          do p=nC+1,nBas-nR
 
            r_2p1h(i,a,b,p) = sqrt(2d0)*NOVV(p,i,b,a) + delta_2p1h(i,a,b,p)*t_2p1h(i,a,b,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2p1h(p,i,a,b) = r_2p1h(p,i,a,b) + 2d0*VOOV(a,k,i,c)*t_2p1h(k,c,b,p)

                do l=nC+1,nO
                  do q=nC+1,nBas-nR

                    r_2p1h(p,i,a,b) = r_2p1h(p,i,a,b) - sqrt(2d0)*t_2p1h(i,a,b,q)*NVOO(q,c,k,l)*t_2h1p(k,l,c,p)

                  end do
                end do

                do d=1,nV-nR
                  do q=nC+1,nBas-nR

                    r_2p1h(p,i,a,b) = r_2p1h(p,i,a,b) - sqrt(2d0)*t_2p1h(i,a,b,q)*NOVV(q,k,d,c)*t_2p1h(k,c,d,p)

                  end do
                end do

              end do
            end do
 
          end do
        end do
      end do
    end do

!   Check convergence 

    Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
  
!   Update amplitudes

    t_2h1p = t_2h1p - r_2h1p/delta_2h1p
    t_2p1h = t_2p1h - r_2p1h/delta_2p1h

!   Compute correlation energy

    eGW(:) = e(:)
    do p=nC+1,nBas-nR

      do i=nC+1,nO
        do j=nC+1,nO
          do a=1,nV-nR

            eGW(p) = eGW(p) + sqrt(2d0)*t_2h1p(i,j,a,p)*NVOO(p,a,i,j)

          end do
        end do
      end do

      do i=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR

            eGW(p) = eGW(p) + sqrt(2d0)*t_2p1h(i,a,b,p)*NOVV(p,i,a,b)

          end do
        end do
      end do

    end do

!   Renormalization factor

    Z(:) = 1d0

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',eGW(nO),'|',eGW(nO+1),'|',Conv,'|'

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

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)'  CCGW calculation  '
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
            '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
  write(*,*)'-------------------------------------------------------------------------------'

  do p=1,nBas
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',e(p)*HaToeV,'|',(eGW(p)-e(p))*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  enddo
  write(*,*)'-------------------------------------------------------------------------------'

!------------------------------------------------------------------------
! EOM section
!------------------------------------------------------------------------

end subroutine CCGW
