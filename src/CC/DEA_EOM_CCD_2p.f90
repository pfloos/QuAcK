subroutine DEA_EOM_CCD_2p(nC,nO,nV,nR,eV,OOVV,VVVV,t)

! DEA-EOM-CCD calculation up to 2p

  implicit none

! Input variables

  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eV(nV)
  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)
  double precision,intent(in)   :: VVVV(nV,nV,nV,nV)
  double precision,intent(in)   :: t(nO,nO,nV,nV)
  
! Local variables

  integer                       :: a,b,c,d,ab,cd
  integer                       :: i,j
  integer                       :: nVV
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: W(:,:,:,:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: Om(:)


! Hello world

  write(*,*)
  write(*,*)'********************'
  write(*,*)'| DEA-EOM-CCD (2p) |'
  write(*,*)'********************'
  write(*,*)

! Size of the EOM Hamiltonian

  nVV = nV*(nV-1)/2

! Memory allocation

  allocate(F(nV,nV),W(nV,nV,nV,nV),H(nVV,nVV),Om(nVV))

! Form one-body terms

  do a=1,nV-nR
    do b=1,nV-nR
 
      F(a,b) = eV(a)*Kronecker_delta(a,b) 

      do i=1,nO-nR
        do j=1,nO-nR
          do c=1,nV-nC
    
            F(a,b) = F(a,b) - 0.5d0*OOVV(i,j,b,c)*t(i,j,a,c)

          end do
        end do
      end do

    end do
  end do

! Form two-body terms

  do a=1,nV-nR
    do b=1,nV-nR
      do c=1,nV-nR
        do d=1,nV-nR
 
          W(a,b,c,d) = VVVV(a,b,c,d)

          do i=1,nO-nC
            do j=i+1,nO-nC
    
              W(a,b,c,d) = W(a,b,c,d) + OOVV(i,j,a,b)*t(i,j,a,b)

            end do
          end do

        end do
      end do
    end do
  end do

! Form EOM Hamiltonian

  ab = 0
  do a=1,nV-nR
   do b=a+1,nV-nR
      ab = ab + 1

      cd = 0
      do c=1,nV-nR
       do d=c+1,nV-nR
          cd = cd + 1

          H(ab,cd) = F(a,c)*Kronecker_delta(b,d) + Kronecker_delta(a,c)*F(b,d) + W(a,b,c,d)

        end do
      end do

    end do
  end do

! Diagonalize EOM Hamiltonian

  if(nVV > 0) call diagonalize_matrix(nVV,H,Om)

! Dump results

  call print_excitation('DEA-EOM-CCD ',3,nVV,Om)

end subroutine 
