subroutine CIS_D(ispin,nBasin,nCin,nOin,nVin,nRin,nSin,maxS,eHF,ERI,w,X)

! Compute the D correction of CIS(D)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBasin
  integer,intent(in)            :: nCin
  integer,intent(in)            :: nOin
  integer,intent(in)            :: nVin
  integer,intent(in)            :: nRin
  integer,intent(in)            :: nSin
  integer,intent(in)            :: maxS
  double precision,intent(in)   :: eHF(nBasin)
  double precision,intent(in)   :: ERI(nBasin,nBasin,nBasin,nBasin)
  double precision,intent(in)   :: w(maxS)
  double precision,intent(in)   :: X(nSin,maxS)

! Local variables

  integer                       :: i,j,k
  integer                       :: a,b,c
  integer                       :: m,ia

  integer                       :: nBas
  integer                       :: nC
  integer                       :: nO
  integer                       :: nV
  integer                       :: nR
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta(:,:,:,:)

  double precision,allocatable  :: OOOV(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVV(:,:,:,:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: X3(:,:)

  double precision,allocatable  :: u(:,:,:,:)
  double precision,allocatable  :: v(:,:)
  double precision,allocatable  :: t(:,:,:,:)
  double precision,allocatable  :: rr(:,:),r(:,:)
  double precision              :: wD

  double precision,external     :: Kronecker_delta

! Spatial to spin orbitals

  nBas = 2*nBasin
  nC   = 2*nCin
  nO   = 2*nOin
  nV   = 2*nVin
  nR   = 2*nRin

  allocate(seHF(nBas),sERI(nBas,nBas,nBas,nBas))

  call spatial_to_spin_MO_energy(nBasin,eHF,nBas,seHF)
  call spatial_to_spin_ERI(nBasin,ERI,nBas,sERI)

! Antysymmetrize ERIs

  allocate(dbERI(nBas,nBas,nBas,nBas))

  call antisymmetrize_ERI(2,nBas,sERI,dbERI)

  deallocate(sERI)

! Form energy denominator

  allocate(eO(nO),eV(nV))
  allocate(delta(nO,nO,nV,nV))

  eO(:) = seHF(1:nO)
  eV(:) = seHF(nO+1:nBas)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta)

  deallocate(seHF,eO,eV)

! Create integral batches

  allocate(OOOV(nO,nO,nO,nV),OOVV(nO,nO,nV,nV),OVVV(nO,nV,nV,nV))

  OOOV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,   1:nO  ,nO+1:nBas)
  OOVV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVVV(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,nO+1:nBas,nO+1:nBas)

  deallocate(dbERI)
 
! Memory allocation

  allocate(t(nO,nO,nV,nV),r(nO,nV),u(nO,nO,nV,nV),v(nO,nV))

  allocate(X1(nV,nV),X2(nO,nO),X3(nO,nV))

! MP2 guess amplitudes

  t(:,:,:,:) = -OOVV(:,:,:,:)/delta(:,:,:,:)

!------------------------------------------------------------------------
! Loop over single excitations
!------------------------------------------------------------------------

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' CIS(D) correction '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20)') '#','CIS (eV)','CIS(D) (eV)','Correction (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do m=1,maxS

    ! Unfold r

    allocate(rr(nOin,nVin))

    ia = 0
    do i=nCin+1,nOin
      do a=1,nVin-nRin
        ia = ia + 1
        rr(i,a) = x(ia,m)
      end do
    end do
 
    if(ispin == 1) then

      do i=nC+1,nO
        do a=1,nV-nR
          r(i,a) = rr((i+1)/2,(a+1)/2)*Kronecker_delta(mod(i,2),mod(a,2))
        end do
      end do

    elseif(ispin == 2) then

      do i=nC+1,nO
        do a=1,nV-nR
          r(i,a) = rr((i+1)/2,(a+1)/2)*Kronecker_delta(mod(i,2),mod(a+1,2))
        end do
      end do

    else 

      print*,'!!! CIS(D) must be for singlet or triplet !!!'
      stop

    end if

    deallocate(rr)

    ! Compute u array

    u(:,:,:,:) = 0d0
 
    do i=nC+1,nO
      do a=1,nV-nR
        do j=nC+1,nO
          do b=1,nV-nR
 
            do c=1,nV-nR
              u(i,j,a,b) = u(i,j,a,b) + OVVV(i,c,a,b)*r(j,c) - OVVV(j,c,a,b)*r(i,c)
            end do
 
            do k=nC+1,nO
              u(i,j,a,b) = u(i,j,a,b) + OOOV(i,j,k,a)*r(k,b) - OOOV(i,j,k,b)*r(k,a)
            end do
 
          end do
        end do
      end do
    end do
  
   ! Compute intermediate arrays

   X1(:,:) = 0d0

  do j=nC+1,nO
    do k=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
          do c=1,nV-nR
            X1(a,b) = X1(a,b) + OOVV(j,k,b,c)*t(j,k,c,a)
          end do
        end do
      end do
    end do
  end do

   X2(:,:) = 0d0

    do i=nC+1,nO
      do j=nC+1,nO
        do k=nC+1,nO
          do b=1,nV-nR
            do c=1,nV-nR
              X2(i,j) = X2(i,j) + OOVV(j,k,b,c)*t(i,k,c,b)
            end do
          end do
        end do
      end do
    end do

   X3(:,:) = 0d0

    do j=nC+1,nO
      do k=nC+1,nO
        do b=1,nV-nR
          do c=1,nV-nR
            X3(k,c) = X3(k,c) + 2d0*OOVV(j,k,b,c)*r(j,b)
          end do
        end do
      end do
    end do

   ! Compute v array

    v(:,:) = 0d0
    
    do i=nC+1,nO
      do a=1,nV-nR

        do b=1,nV-nR
          v(i,a) = v(i,a) + r(i,b)*X1(a,b) 
        end do

        do j=nC+1,nO
          v(i,a) = v(i,a) + r(j,a)*X2(i,j)
        end do

        do k=nC+1,nO
          do c=1,nV-nR
            v(i,a) = v(i,a) + X3(k,c)*t(i,k,a,c)
          end do
        end do

      end do
    end do

    v(:,:) = 0.5d0*v(:,:)

    ! Compute CIS(D) correction to CIS excitation energies

    wD = 0d0
 
    do i=nC+1,nO
      do a=1,nV-nR
        do j=nC+1,nO
          do b=1,nV-nR

            wD = wD - 0.25d0*u(i,j,a,b)**2/(delta(i,j,a,b) - w(m)) 

          end do
        end do
        wD = wD + r(i,a)*v(i,a)
      end do
    end do
    wD = 0.5d0*wD

    ! Flush results
 
    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6)') m,w(m)*HaToeV,(w(m)+wD)*HaToeV,wD*HaToeV

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*)
!------------------------------------------------------------------------
! End of loop over single excitations
!------------------------------------------------------------------------

end subroutine 
