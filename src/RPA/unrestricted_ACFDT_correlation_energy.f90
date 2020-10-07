subroutine unrestricted_ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt, & 
                                                 ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,XmY,EcAC)

! Compute the correlation energy via the adiabatic connection formula

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  logical,intent(in)            :: exchange_kernel
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: ia,jb
  double precision              :: delta_Kx
  double precision,allocatable  :: Ap(:,:)
  double precision,allocatable  :: Bp(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: EcAC

! Exchange kernel

  delta_Kx = 0d0
  if(exchange_kernel) delta_Kx = 1d0
  
! Memory allocation

  allocate(Ap(nSt,nSt),Bp(nSt,nSt),X(nSt,nSt),Y(nSt,nSt))

! Compute Aiajb = (ia|bj) and Biajb = (ia|jb)

! Initialization

  Ap(:,:) = 0d0
  Bp(:,:) = 0d0

!-----------------------------------------------
! Build kernel for spin-conserving transitions
!-----------------------------------------------

  if(ispin == 1) then

    ! aaaa block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            Ap(ia,jb) = ERI_aaaa(i,b,a,j) - delta_Kx*ERI_aaaa(i,b,j,a)
            Bp(ia,jb) = ERI_aaaa(i,j,a,b) - delta_Kx*ERI_aaaa(i,j,b,a)

          end  do
        end  do
      end  do
    end  do

    ! aabb block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            Ap(ia,nSa+jb) = ERI_aabb(i,b,a,j)
            Bp(ia,nSa+jb) = ERI_aabb(i,j,a,b)


          end  do
        end  do
      end  do
    end  do

    ! bbaa block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            Ap(nSa+ia,jb) = ERI_aabb(b,i,j,a)
            Bp(nSa+ia,jb) = ERI_aabb(j,i,b,a)

          end  do
        end  do
      end  do
    end  do

    ! bbbb block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            Ap(nSa+ia,nSa+jb) = ERI_bbbb(i,b,a,j) - delta_Kx*ERI_bbbb(i,b,j,a)
            Bp(nSa+ia,nSa+jb) = ERI_bbbb(i,j,a,b) - delta_Kx*ERI_bbbb(i,j,b,a)


          end  do
        end  do
      end  do
    end  do

  end if

!-----------------------------------------------
! Build A matrix for spin-flip transitions
!-----------------------------------------------

  if(ispin == 2) then

    ! abab block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1
            Ap(ia,jb) = - delta_Kx*ERI_aabb(i,b,j,a)

          end  do
        end  do
      end  do
    end  do

    ! baba block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            Ap(nSa+ia,nSa+jb) = - delta_Kx*ERI_aabb(b,i,a,j)

          end  do
        end  do
      end  do
    end  do

    ! abba block

    ia = 0
    do i=nC(1)+1,nO(1)
      do a=nO(2)+1,nBas-nR(2)
        ia = ia + 1
        jb = 0
        do j=nC(2)+1,nO(2)
          do b=nO(1)+1,nBas-nR(1)
            jb = jb + 1

            Bp(ia,nSa+jb) = - delta_Kx*ERI_aabb(i,j,b,a)

          end  do
        end  do
      end  do
    end  do

    ! baab block

    ia = 0
    do i=nC(2)+1,nO(2)
      do a=nO(1)+1,nBas-nR(1)
        ia = ia + 1
        jb = 0
        do j=nC(1)+1,nO(1)
          do b=nO(2)+1,nBas-nR(2)
            jb = jb + 1

            Bp(nSa+ia,jb) = - delta_Kx*ERI_aabb(j,i,a,b)

          end  do
        end  do
      end  do
    end  do

  end if

! Compute Tr(K x P_lambda)

  X(:,:) = 0.5d0*(XpY(:,:) + XmY(:,:))
  Y(:,:) = 0.5d0*(XpY(:,:) - XmY(:,:))

  EcAC = trace_matrix(nSt,matmul(X,matmul(Bp,transpose(Y))) + matmul(Y,matmul(Bp,transpose(X)))) &
       + trace_matrix(nSt,matmul(X,matmul(Ap,transpose(X))) + matmul(Y,matmul(Ap,transpose(Y)))) &
       - trace_matrix(nSt,Ap)

end subroutine unrestricted_ACFDT_correlation_energy
