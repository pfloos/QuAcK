subroutine ccRG0W0(maxSCF,thresh,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! CC-based GW module

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer                       :: p,q
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  integer                       :: nSCF
  double precision              :: Conv

  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VOOV(:,:,:,:)

  double precision,allocatable  :: delta_2h1p(:,:,:)
  double precision,allocatable  :: delta_2p1h(:,:,:)

  double precision,allocatable  :: V_2h1p(:,:,:)
  double precision,allocatable  :: V_2p1h(:,:,:)

  double precision,allocatable  :: r_2h1p(:,:,:)
  double precision,allocatable  :: r_2p1h(:,:,:)

  double precision,allocatable  :: t_2h1p(:,:,:)
  double precision,allocatable  :: t_2p1h(:,:,:)

  double precision              :: x_2h1p
  double precision              :: x_2p1h

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)


! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* CC-based G0W0 Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Create integral batches

  allocate(OVVO(nO,nV,nV,nO),VOOV(nV,nO,nO,nV))

  OVVO(:,:,:,:) = ERI(   1:nO   ,nO+1:nOrb,nO+1:nOrb,   1:nO  )
  VOOV(:,:,:,:) = ERI(nO+1:nOrb ,   1:nO  ,   1:nO  ,nO+1:nOrb)
 
! Form energy denominator and guess amplitudes

  allocate(delta_2h1p(nO,nO,nV),delta_2p1h(nO,nV,nV))
  allocate(V_2h1p(nO,nO,nV),V_2p1h(nO,nV,nV))
  allocate(t_2h1p(nO,nO,nV),t_2p1h(nO,nV,nV))
  allocate(r_2h1p(nO,nO,nV),r_2p1h(nO,nV,nV))
  allocate(eGW(nOrb),Z(nOrb))

! Initialization

  eGW(:) = eHF(:)

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

    ! Initialization
 
    Conv = 1d0
    nSCF =  0
 
    t_2h1p(:,:,:) = 0d0
    t_2p1h(:,:,:) = 0d0
 
    ! Compute energy differences
  
    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR
  
          delta_2h1p(i,j,a) = eHF(i) + eHF(j) - eHF(nO+a) - eHF(p)
  
        end do
      end do
    end do
  
    do i=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
  
          delta_2p1h(i,a,b) = eHF(nO+a) + eHF(nO+b) - eHF(i) - eHF(p)
  
        end do
      end do
    end do

    ! Compute V2h1p and V2p1h

    do k=nC+1,nO
      do l=nC+1,nO
        do c=1,nV-nR
 
          V_2h1p(k,l,c) = sqrt(2d0)*ERI(p,nO+c,k,l)
 
        end do
      end do
    end do
 
    do k=nC+1,nO
      do c=1,nV-nR
        do d=1,nV-nR
 
          V_2p1h(k,c,d) = sqrt(2d0)*ERI(p,k,nO+d,nO+c)
 
        end do
      end do
    end do

   !----------------------!
   ! Loop over amplitudes !
   !----------------------!

    write(*,*)
    write(*,*)'----------------------------------------------'
    write(*,*)'| CC-based G0W0 calculation                  |'
    write(*,*)'----------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
              '|','#','|','HF','|','G0W0','|','Conv','|'
    write(*,*)'----------------------------------------------'
 
    do while(Conv > thresh .and. nSCF < maxSCF)
 
     ! Increment 
 
      nSCF = nSCF + 1
 
      !  Compute intermediates
 
      x_2h1p = 0d0 
 
      do k=nC+1,nO
        do l=nC+1,nO
          do c=1,nV-nR
      
            x_2h1p = x_2h1p + V_2h1p(k,l,c)*t_2h1p(k,l,c)
      
          end do
        end do
      end do
     
      x_2p1h = 0d0 
 
      do k=nC+1,nO
        do c=1,nV-nR
          do d=1,nV-nR
      
            x_2p1h = x_2p1h + V_2p1h(k,c,d)*t_2p1h(k,c,d)
     
          end do
        end do
      end do
     
      ! Compute residual for 2h1p sector
 
      do i=nC+1,nO
        do j=nC+1,nO
          do a=1,nV-nR
 
            r_2h1p(i,j,a) = V_2h1p(i,j,a) + delta_2h1p(i,j,a)*t_2h1p(i,j,a)
 
            do k=nC+1,nO
              do c=1,nV-nR
 
                r_2h1p(i,j,a) = r_2h1p(i,j,a) - 2d0*OVVO(j,c,a,k)*t_2h1p(i,k,c)
 
              end do
            end do

            r_2h1p(i,j,a) = r_2h1p(i,j,a) - t_2h1p(i,j,a)*x_2h1p - t_2h1p(i,j,a)*x_2p1h
 
          end do
        end do
      end do
 
      ! Compute residual for 2p1h sector
 
      do i=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR
 
            r_2p1h(i,a,b) = V_2p1h(i,a,b) + delta_2p1h(i,a,b)*t_2p1h(i,a,b)
 
            do k=nC+1,nO
              do c=1,nV-nR
 
                r_2p1h(i,a,b) = r_2p1h(i,a,b) + 2d0*VOOV(a,k,i,c)*t_2p1h(k,c,b)
 
              end do
            end do
 
            r_2p1h(i,a,b) = r_2p1h(i,a,b) - t_2p1h(i,a,b)*x_2h1p - t_2p1h(i,a,b)*x_2p1h
 
          end do
        end do
      end do
  
      !  Check convergence 
 
      Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
    
      ! Update amplitudes
 
      t_2h1p(:,:,:) = t_2h1p(:,:,:) - r_2h1p(:,:,:)/delta_2h1p(:,:,:)
      t_2p1h(:,:,:) = t_2p1h(:,:,:) - r_2p1h(:,:,:)/delta_2p1h(:,:,:)
 
      ! Compute self-energy
 
      eGW(p) = eHF(p)
 
      do i=nC+1,nO
        do j=nC+1,nO
          do a=1,nV-nR
 
            eGW(p) = eGW(p) + V_2h1p(i,j,a)*t_2h1p(i,j,a)
  
          end do
        end do
      end do
 
      do i=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR
 
            eGW(p) = eGW(p) + V_2p1h(i,a,b)*t_2p1h(i,a,b)
  
          end do
        end do
      end do
 
      ! Renormalization factor
 
      Z(:) = 1d0
 
      ! Dump results
 
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
        '|',nSCF,'|',eHF(p)*HaToeV,'|',eGW(p)*HaToeV,'|',Conv,'|'
 
    end do

    write(*,*)'----------------------------------------------'
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
 
    end if

    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)' CC-G0W0 calculation '
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
              '|','#','|','e_HF (eV)','|','Sig_c (eV)','|','Z','|','e_QP (eV)','|'
    write(*,*)'-------------------------------------------------------------------------------'
 
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',(eGW(p)-eHF(p))*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
    write(*,*)'-------------------------------------------------------------------------------'

  end do

end subroutine 
