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

  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VOOV(:,:,:,:)

  double precision,allocatable  :: delta_2h1p(:,:,:,:)
  double precision,allocatable  :: delta_2p1h(:,:,:,:)

  double precision,allocatable  :: V_2h1p(:,:,:,:)
  double precision,allocatable  :: V_2p1h(:,:,:,:)

  double precision,allocatable  :: r_2h1p(:,:,:,:)
  double precision,allocatable  :: r_2p1h(:,:,:,:)

  double precision,allocatable  :: t_2h1p(:,:,:,:)
  double precision,allocatable  :: t_2p1h(:,:,:,:)

  double precision,allocatable  :: x_2h1p(:,:)
  double precision,allocatable  :: x_2p1h(:,:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: SigGW(:,:)
  double precision,allocatable  :: cGW(:,:)
  double precision,allocatable  :: Z(:)

  integer,allocatable           :: order(:)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'|     CCGW calculation      |'
  write(*,*)'*****************************'
  write(*,*)

! Create integral batches

  allocate(OVVO(nO,nV,nV,nO),VOOV(nV,nO,nO,nV))

  OVVO(:,:,:,:) = ERI(   1:nO   ,nO+1:nBas,nO+1:nBas,   1:nO  )
  VOOV(:,:,:,:) = ERI(nO+1:nBas ,   1:nO  ,   1:nO  ,nO+1:nBas)
 
! Form energy denominator and guess amplitudes

  allocate(eO(nO),eV(nV))
  allocate(delta_2h1p(nO,nO,nV,nBas),delta_2p1h(nO,nV,nV,nBas))
  allocate(V_2h1p(nBas,nO,nO,nV),V_2p1h(nBas,nO,nV,nV))
  allocate(t_2h1p(nO,nO,nV,nBas),t_2p1h(nO,nV,nV,nBas))
  allocate(x_2h1p(nBas,nBas),x_2p1h(nBas,nBas))

  eO(:) = e(1:nO)
  eV(:) = e(nO+1:nBas)

  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do p=nC+1,nBas-nR

          delta_2h1p(i,j,a,p) = eO(i) + eO(j) - eV(a) - e(p)
          V_2h1p(p,i,j,a) = sqrt(2d0)*ERI(p,nO+a,j,i)

        end do
      end do
    end do
  end do

  do i=nC+1,nO
    do a=1,nV-nR
      do b=1,nV-nR
        do p=nC+1,nBas-nR

          delta_2p1h(i,a,b,p) = eV(a) + eV(b) - eO(i) - e(p)
          V_2p1h(p,i,a,b) = sqrt(2d0)*ERI(p,i,nO+b,nO+a)

        end do
      end do
    end do
  end do

! Initialization

  allocate(r_2h1p(nO,nO,nV,nBas),r_2p1h(nO,nV,nV,nBas))
  allocate(eGW(nBas),SigGW(nBas,nBas),cGW(nBas,nBas),Z(nBas))
  allocate(order(nBas))

  Conv = 1d0
  nSCF = 0

  t_2h1p(:,:,:,:) = 0d0
  t_2p1h(:,:,:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------'
  write(*,*)'| CCGW calculation                           |'
  write(*,*)'----------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','HOMO','|','LUMO','|','Conv','|'
  write(*,*)'----------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Compute intermediates

    x_2h1p(:,:) = 0d0 

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR

        do k=nC+1,nO
          do l=nC+1,nO
            do c=1,nV-nR
      
              x_2h1p(p,q) = x_2h1p(p,q) + V_2h1p(q,k,l,c)*t_2h1p(k,l,c,p)
      
            end do
          end do
        end do

      end do
    end do
   
    x_2p1h(:,:) = 0d0 

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR

        do k=nC+1,nO
          do c=1,nV-nR
            do d=1,nV-nR
      
              x_2p1h(p,q) = x_2p1h(p,q) + V_2p1h(q,k,c,d)*t_2p1h(k,c,d,p)
     
          end do
        end do
      end do

      end do
    end do
   
!   Compute residual for 2h1p sector

    do i=nC+1,nO
      do j=nC+1,nO
        do a=1,nV-nR

          do p=nC+1,nBas-nR

            r_2h1p(i,j,a,p) = V_2h1p(p,i,j,a) + delta_2h1p(i,j,a,p)*t_2h1p(i,j,a,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - 2d0*OVVO(j,c,a,k)*t_2h1p(i,k,c,p)

              end do
            end do

            do q=nC+1,nBas-nR

              r_2h1p(i,j,a,p) = r_2h1p(i,j,a,p) - t_2h1p(i,j,a,q)*x_2h1p(p,q) - t_2h1p(i,j,a,q)*x_2p1h(p,q)

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

            r_2p1h(i,a,b,p) = V_2p1h(p,i,a,b) + delta_2p1h(i,a,b,p)*t_2p1h(i,a,b,p)

            do k=nC+1,nO
              do c=1,nV-nR

                r_2p1h(i,a,b,p) = r_2p1h(i,a,b,p) + 2d0*VOOV(a,k,i,c)*t_2p1h(k,c,b,p)

              end do
            end do

            do q=nC+1,nBas-nR

              r_2p1h(i,a,b,p) = r_2p1h(i,a,b,p) - t_2p1h(i,a,b,q)*x_2h1p(p,q) - t_2p1h(i,a,b,q)*x_2p1h(p,q)

            end do
 
          end do

        end do
      end do
    end do
 
!   Check convergence 

    Conv = max(maxval(abs(r_2h1p)),maxval(abs(r_2p1h)))
  
!   Update amplitudes

    t_2h1p(:,:,:,:) = t_2h1p(:,:,:,:) - r_2h1p(:,:,:,:)/delta_2h1p(:,:,:,:)
    t_2p1h(:,:,:,:) = t_2p1h(:,:,:,:) - r_2p1h(:,:,:,:)/delta_2p1h(:,:,:,:)

!   Compute correlation energy

    SigGW(:,:) = 0d0

    do p=nC+1,nBas-nR

      SigGW(p,p) = SigGW(p,p) + e(p)

      do q=nC+1,nBas-nR

        do i=nC+1,nO
          do j=nC+1,nO
            do a=1,nV-nR

              SigGW(p,q) = SigGW(p,q) + V_2h1p(p,i,j,a)*t_2h1p(i,j,a,q)
 
            end do
          end do
        end do

        do i=nC+1,nO
          do a=1,nV-nR
            do b=1,nV-nR

              SigGW(p,q) = SigGW(p,q) + V_2p1h(p,i,a,b)*t_2p1h(i,a,b,q)
 
            end do
          end do
        end do

      end do
    end do

    call diagonalize_general_matrix(nBas,SigGW,eGW,cGW)

    do p=1,nBas
      order(p) = p
    end do

    call quick_sort(eGW,order,nBas)
    call set_order(cGW,order,nBas,nBas)

!   Renormalization factor

    Z(:) = 1d0

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',eGW(nO)*HaToeV,'|',eGW(nO+1)*HaToeV,'|',Conv,'|'

  enddo
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

end subroutine CCGW
