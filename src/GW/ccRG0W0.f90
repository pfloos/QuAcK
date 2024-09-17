subroutine ccRG0W0(maxSCF,thresh,max_diis,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

! CC-based GW module

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis

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

  integer                       :: p,q,r,s
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

  integer                       :: nSCF
  double precision              :: Conv

  double precision              :: x

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)

  double precision,allocatable  :: del(:,:,:)
  double precision,allocatable  :: res(:,:,:)
  double precision,allocatable  :: amp(:,:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: r_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)


! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* CC-based G0W0 Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Form energy denominator and guess amplitudes

  allocate(del(nOrb,nOrb,nOrb))
  allocate(res(nOrb,nOrb,nOrb))
  allocate(amp(nOrb,nOrb,nOrb))

  allocate(eGW(nOrb),Z(nOrb))

 allocate(r_diis(nOrb**3,max_diis))
 allocate(t_diis(nOrb**3,max_diis))

! Initialization

  eGW(:) = eHF(:)
  Z(:) = 1d0

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  do p=nO,nO

    ! Initialization
 
    Conv = 1d0
    nSCF =  0

    n_diis      = 0
    t_diis(:,:) = 0d0
    r_diis(:,:) = 0d0
    rcond       = 0d0

    amp(:,:,:) = 0d0
    res(:,:,:) = 0d0
    del(:,:,:) = huge(1d0)
 
    ! Compute energy differences
 
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR
  
          del(i,j,a) = eHF(i) + eHF(j) - eHF(a) - eHF(p)
  
        end do
      end do
    end do
  
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
  
          del(b,a,i) = eHF(a) + eHF(b) - eHF(i) - eHF(p)
  
        end do
      end do
    end do

   !----------------------!
   ! Loop over amplitudes !
   !----------------------!

    write(*,*)
    write(*,*)'-------------------------------------------------------------'
    write(*,*)'| CC-based G0W0 calculation                                 |'
    write(*,*)'-------------------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X,A15,1X,A1,1X)') &
              '|','#','|','HF','|','G0W0','|','Conv','|'
    write(*,*)'-------------------------------------------------------------'
 
    do while(Conv > thresh .and. nSCF < maxSCF)
 
     ! Increment 
 
      nSCF = nSCF + 1
 
      !  Compute intermediates x_2h1p and x_2p1h
 
      x = 0d0 
 
      do q=nC+1,nOrb-nR
        do r=nC+1,nOrb-nR
          do s=nC+1,nOrb-nR

            x = x + sqrt(2d0)*ERI(p,s,q,r)*amp(q,r,s)
     
          end do
        end do
      end do

      ! Compute residual for 2h1p sector
 
      do i=nC+1,nO
        do j=nC+1,nO
          do a=nO+1,nOrb-nR
 
            res(i,j,a) = sqrt(2d0)*ERI(p,a,i,j) + (del(i,j,a) - x)*amp(i,j,a)
 
            do k=nC+1,nO
              do c=nO+1,nOrb-nR
 
                res(i,j,a) = res(i,j,a) - 2d0*ERI(j,c,a,k)*amp(i,k,c)
 
              end do
            end do
 
          end do
        end do
      end do
 
      ! Compute residual for 2p1h sector
 
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          do b=nO+1,nOrb-nR
 
            res(b,a,i) = sqrt(2d0)*ERI(p,i,b,a) + (del(b,a,i) - x)*amp(b,a,i)
 
            do k=nC+1,nO
              do c=nO+1,nOrb-nR
 
                res(b,a,i) = res(b,a,i) + 2d0*ERI(a,k,i,c)*amp(b,c,k) 
 
              end do
            end do
 
          end do
        end do
      end do
  
      !  Check convergence 
 
      Conv = maxval(abs(res))
    
      ! Update amplitudes

      amp(:,:,:) = amp(:,:,:) - res(:,:,:)/del(:,:,:)
 
      ! DIIS extrapolation

      if(max_diis > 1) then
     
        n_diis = min(n_diis+1,max_diis)
        call DIIS_extrapolation(rcond,nOrb**3,nOrb**3,n_diis,r_diis,t_diis,res,amp)
     
      end if

      ! Compute quasiparticle energy
 
      eGW(p) = eHF(p)

      do q=nC+1,nOrb-nR
        do r=nC+1,nOrb-nR
          do s=nC+1,nOrb-nR
 
            eGW(p) = eGW(p) + sqrt(2d0)*ERI(p,s,q,r)*amp(q,r,s)
  
          end do
        end do
      end do

      ! Dump results
 
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X)') &
        '|',nSCF,'|',eHF(p)*HaToeV,'|',eGW(p)*HaToeV,'|',Conv,'|'
 
    end do

    write(*,*)'-------------------------------------------------------------'
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
 
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X,F15.10,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',(eGW(p)-eHF(p))*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
    write(*,*)'-------------------------------------------------------------------------------'

  end do

end subroutine 
