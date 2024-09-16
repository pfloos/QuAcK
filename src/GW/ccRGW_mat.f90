subroutine ccRGW_mat(maxSCF,thresh,nBas,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF)

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

  double precision,allocatable  :: delta_2h1p(:,:)
  double precision,allocatable  :: delta_2p1h(:,:)

  double precision,allocatable  :: C_2h1p(:,:)
  double precision,allocatable  :: C_2p1h(:,:)

  double precision,allocatable  :: V_2h1p(:,:)
  double precision,allocatable  :: V_2p1h(:,:)

  double precision,allocatable  :: r_2h1p(:,:)
  double precision,allocatable  :: r_2p1h(:,:)

  double precision,allocatable  :: t_2h1p(:,:)
  double precision,allocatable  :: t_2p1h(:,:)

  double precision,allocatable  :: SigGW(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: F(:,:)
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

  allocate(delta_2h1p(n2h1p,nOrb),delta_2p1h(n2p1h,nOrb))
  allocate(C_2h1p(n2h1p,n2h1p),C_2p1h(n2p1h,n2p1h))
  allocate(V_2h1p(nOrb,n2h1p),V_2p1h(nOrb,n2p1h))
  allocate(t_2h1p(n2h1p,nOrb),t_2p1h(n2p1h,nOrb))
  allocate(r_2h1p(n2h1p,nOrb),r_2p1h(n2p1h,nOrb))
  allocate(F(nOrb,nOrb),eGW(nOrb),SigGW(nOrb,nOrb),Z(nOrb))

  F(:,:) = 0d0
  do p=nC+1,nOrb-nR
    F(p,p) = eHF(p) 
  end do

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

  ! Initialization

  Conv = 1d0
  nSCF =  0

  t_2h1p(:,:) = 0d0
  t_2p1h(:,:) = 0d0

  ! Compute energy differences

  ija = 0
  do i=nC+1,nO
    do j=nC+1,nO
      do a=nO+1,nOrb-nR
        ija = ija + 1
  
        do p=nC+1,nOrb-nR

          delta_2h1p(ija,p) = eHF(i) + eHF(j) - eHF(a) - eGW(p)

        end do

      end do
    end do
  end do

  iab = 0
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=nO+1,nOrb-nR
        iab = iab + 1

        do p=nC+1,nOrb-nR

          delta_2p1h(iab,p) = eHF(a) + eHF(b) - eHF(i) - eGW(p)

        end do

      end do
    end do
  end do

  klc = 0
  do k=nC+1,nO
    do l=nC+1,nO
      do c=nO+1,nOrb-nR
        klc = klc + 1

        do p=nC+1,nOrb-nR

          V_2h1p(p,klc) = sqrt(2d0)*ERI(p,c,k,l)

        end do

      end do
    end do
  end do

  kcd = 0
  do k=nC+1,nO
    do c=nO+1,nOrb-nR
      do d=nO+1,nOrb-nR
        kcd = kcd + 1

        do p=nC+1,nOrb-nR

          V_2p1h(p,kcd) = sqrt(2d0)*ERI(p,k,d,c)

        end do

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

    r_2h1p = transpose(V_2h1p) + matmul(C_2h1p,t_2h1p) - matmul(t_2h1p,F) &
           - matmul(matmul(t_2h1p,V_2h1p),t_2h1p) - matmul(matmul(t_2h1p,V_2p1h),t_2p1h)

    ! Compute residual for 2p1h sector

    r_2p1h = transpose(V_2p1h) + matmul(C_2p1h,t_2p1h) - matmul(t_2p1h,F) &
           - matmul(matmul(t_2p1h,V_2h1p),t_2h1p) - matmul(matmul(t_2p1h,V_2p1h),t_2p1h)

    !  Check convergence 

    Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
  
    ! Update amplitudes

    t_2h1p(:,:) = t_2h1p(:,:) - r_2h1p(:,:)/delta_2h1p(:,:)
    t_2p1h(:,:) = t_2p1h(:,:) - r_2p1h(:,:)/delta_2p1h(:,:)

    ! Compute quasiparticle energies

    SigGW(:,:) = F(:,:) + matmul(V_2h1p,t_2h1p) + matmul(V_2p1h,t_2p1h)

    call diagonalize_matrix(nOrb,SigGW,eGW)

    ! Renormalization factor

    Z(:) = 0d0

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
    '|',nSCF,'|',eGW(nO)*HaToeV,'|',eGW(nO+1)*HaToeV,'|',Conv,'|'

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

  do p=nC+1,nOrb-nR
    write(*,'(1X,A1,1X,I3,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X,F15.6,1X,A1,1X)') &
    '|',p,'|',eHF(p)*HaToeV,'|',(eGW(p)-eHF(p))*HaToeV,'|',Z(p),'|',eGW(p)*HaToeV,'|'
  end do
  write(*,*)'-------------------------------------------------------------------------------'

end subroutine 
