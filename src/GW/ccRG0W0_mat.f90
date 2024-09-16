subroutine ccRG0W0_mat(maxSCF,thresh,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

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

  double precision,allocatable  :: delta_2h1p(:)
  double precision,allocatable  :: delta_2p1h(:)

  double precision,allocatable  :: C_2h1p(:,:)
  double precision,allocatable  :: C_2p1h(:,:)

  double precision,allocatable  :: V_2h1p(:)
  double precision,allocatable  :: V_2p1h(:)

  double precision,allocatable  :: r_2h1p(:)
  double precision,allocatable  :: r_2p1h(:)

  double precision,allocatable  :: t_2h1p(:)
  double precision,allocatable  :: t_2p1h(:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: Z(:)

  double precision,external     :: Kronecker_delta
  
  integer                       :: n1h,n1p,n1h1p,n2h1p,n2p1h
  integer                       :: ija,klc,iab,kcd


! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* CC-based G0W0 Calculation *'
  write(*,*)'*****************************'
  write(*,*)
 
! Form energy denominator and guess amplitudes

  n1h = nO
  n1p = nV
  n1h1p = n1h*n1p
  n2h1p = n1h*n1h1p
  n2p1h = n1p*n1h1p

  allocate(delta_2h1p(n2h1p),delta_2p1h(n2p1h))
  allocate(C_2h1p(n2h1p,n2h1p),C_2p1h(n2p1h,n2p1h))
  allocate(V_2h1p(n2h1p),V_2p1h(n2p1h))
  allocate(t_2h1p(n2h1p),t_2p1h(n2p1h))
  allocate(r_2h1p(n2h1p),r_2p1h(n2p1h))
  allocate(eGW(nOrb),Z(nOrb))

  ! Compute C2h1p and C2p1h

  ija = 0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1

        klc = 0
        do k=nC+1,nO
          do l=nC+1,nO
            do c=nO+1,nOrb-nR
              klc = klc + 1

              C_2h1p(ija,klc) = ((eHF(i) + eHF(j) - eHF(a))*Kronecker_delta(j,l)*Kronecker_delta(a,c) &
                              - 2d0*ERI(j,c,a,l))*Kronecker_delta(i,k)

            end do
          end do
        end do

      end do
    end do
  end do

  iab = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=nO+1,nOrb-nR
        iab = iab + 1

        kcd = 0
        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=nO+1,nOrb-nR
              kcd = kcd + 1

              C_2p1h(iab,kcd) = ((eHF(a) + eHF(b) - eHF(i))*Kronecker_delta(i,k)*Kronecker_delta(a,c) &
                              + 2d0*ERI(a,k,i,c))*Kronecker_delta(b,d)

            end do
          end do
        end do

      end do
    end do
  end do

!-------------------------!
! Main loop over orbitals !
!-------------------------!

  eGW(:) = eHF(:)

  do p=nO,nO

    ! Initialization
 
    Conv = 1d0
    nSCF =  0
 
    t_2h1p(:) = 0d0
    t_2p1h(:) = 0d0
 
    ! Compute energy differences
 
    ija = 0
    do i=nC+1,nO
      do j=nC+1,nO
        do a=nO+1,nOrb-nR
          ija = ija + 1
 
          delta_2h1p(ija) = eHF(i) + eHF(j) - eHF(a) - eHF(p)
 
        end do
      end do
    end do
 
    iab = 0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        do b=nO+1,nOrb-nR
          iab = iab + 1
 
          delta_2p1h(iab) = eHF(a) + eHF(b) - eHF(i) - eHF(p)
 
        end do
      end do
    end do

    klc = 0
    do k=nC+1,nO
      do l=nC+1,nO
        do c=nO+1,nOrb-nR
          klc = klc + 1
 
          V_2h1p(klc) = sqrt(2d0)*ERI(p,c,k,l)
 
        end do
      end do
    end do
 
    kcd = 0
    do k=nC+1,nO
      do c=nO+1,nOrb-nR
        do d=nO+1,nOrb-nR
          kcd = kcd + 1
 
          V_2p1h(kcd) = sqrt(2d0)*ERI(p,k,d,c)
 
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

      ! Compute residual for 2h1p sector
 
      r_2h1p = V_2h1p + matmul(C_2h1p,t_2h1p) - t_2h1p*eHF(p) &
             - dot_product(t_2h1p,V_2h1p)*t_2h1p - t_2h1p*dot_product(V_2p1h,t_2p1h)

      ! Compute residual for 2p1h sector
 
      r_2p1h = V_2p1h + matmul(C_2p1h,t_2p1h) - t_2p1h*eHF(p) &
             - t_2p1h*dot_product(V_2h1p,t_2h1p) - dot_product(t_2p1h,V_2p1h)*t_2p1h
  
      !  Check convergence 
 
      Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
    
      ! Update amplitudes
 
      t_2h1p(:) = t_2h1p(:) - r_2h1p(:)/delta_2h1p(:)
      t_2p1h(:) = t_2p1h(:) - r_2p1h(:)/delta_2p1h(:)
 
      ! Compute quasiparticle energies
 
      eGW(p) = eHF(p) + dot_product(V_2h1p,t_2h1p) + dot_product(V_2p1h,t_2p1h)
 
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
