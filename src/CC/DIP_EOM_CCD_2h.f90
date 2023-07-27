subroutine DIP_EOM_CCD_2h(nC,nO,nV,nR,eO,OOVV,OOOO,t)

! DIP-EOM-CCD calculation up to 2h

  implicit none

! Input variables

  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eO(nO)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: OOOO(nO,nO,nO,nO)
  double precision,intent(in)   :: t(nO,nO,nV,nV)
  
! Local variables

  integer                       :: i,j,k,l,ij,kl
  integer                       :: a,b
  integer                       :: nOO
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: W(:,:,:,:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: Om(:)


! Hello world

  write(*,*)
  write(*,*)'********************'
  write(*,*)'| DIP-EOM-CCD (2h) |'
  write(*,*)'********************'
  write(*,*)

! Size of the EOM Hamiltonian

  nOO = nO*(nO-1)/2

! Memory allocation

  allocate(F(nO,nO),W(nO,nO,nO,nO),H(nOO,nOO),Om(nOO))

! Form one-body terms

  do i=1,nO-nC
    do j=1,nO-nC
 
      F(i,j) = eO(i)*Kronecker_delta(i,j) 

      do k=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
    
            F(i,j) = F(i,j) + 0.5d0*OOVV(i,k,a,b)*t(j,k,a,b)

          end do
        end do
      end do

    end do
  end do

! Form two-body terms

  do i=1,nO-nC
    do j=1,nO-nC
      do k=1,nO-nC
        do l=1,nO-nC
 
          W(i,j,k,l) = OOOO(i,j,k,l)

          do a=1,nV-nR
            do b=a+1,nV-nR
    
              W(i,j,k,l) = W(i,j,k,l) + OOVV(i,j,a,b)*t(i,j,a,b)

            end do
          end do

        end do
      end do
    end do
  end do

! Form EOM Hamiltonian

  ij = 0
  do i=1,nO-nC
   do j=i+1,nO-nC
      ij = ij + 1

      kl = 0
      do k=1,nO-nC
       do l=k+1,nO-nC
          kl = kl + 1

          H(ij,kl) = - F(i,k)*Kronecker_delta(j,l) - Kronecker_delta(i,k)*F(j,l) + W(i,j,k,l)

        end do
      end do

    end do
  end do

! Diagonalize EOM Hamiltonian

  if(nOO > 0) call diagonalize_matrix(nOO,H,Om)

! Dump results

  call print_excitation('DIP-EOM-CCD ',3,nOO,Om)

end subroutine 
