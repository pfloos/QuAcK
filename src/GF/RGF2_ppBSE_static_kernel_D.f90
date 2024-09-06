subroutine RGF2_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOO,lambda,ERI,eGF,KD_sta)

! Compute the resonant part of the static BSE@GF2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: i,j,k,l,m
  integer                       :: e
  integer                       :: ij,kl

! Output variables

  double precision,intent(out)  :: KD_sta(nOO,nOO)

! Initialization

  KD_sta(:,:) = 0d0

! Second-order correlation kernel for the block D of the singlet manifold

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
      do j=i,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = - eGF(e) + eGF(m)
                num = 2d0*ERI(i,e,k,m)*ERI(j,m,l,e) -     ERI(i,e,k,m)*ERI(j,m,e,l) & 
                    -     ERI(i,e,m,k)*ERI(j,m,l,e) -     ERI(i,e,m,k)*ERI(j,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)
            
                dem = - eGF(e) + eGF(m)
                num = 2d0*ERI(j,e,k,m)*ERI(i,m,l,e) -     ERI(j,e,k,m)*ERI(i,m,e,l) & 
                    -     ERI(j,e,m,k)*ERI(i,m,l,e) -     ERI(j,e,m,k)*ERI(i,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)
            
              end do
            end do

            KD_sta(ij,kl) = 2d0*lambda*KD_sta(ij,kl)/sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block D of the triplet manifold

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = - eGF(e) + eGF(m)
                num = 2d0*ERI(i,e,k,m)*ERI(j,m,l,e) -     ERI(i,e,k,m)*ERI(j,m,e,l) &
                    -     ERI(i,e,m,k)*ERI(j,m,l,e) +     ERI(i,e,m,k)*ERI(j,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) + 2d0*num*dem/(dem**2 + eta**2)

                dem = - eGF(e) + eGF(m)
                num = 2d0*ERI(j,e,k,m)*ERI(i,m,l,e) -     ERI(j,e,k,m)*ERI(i,m,e,l) &
                    -     ERI(j,e,m,k)*ERI(i,m,l,e) +     ERI(j,e,m,k)*ERI(i,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) - 2d0*num*dem/(dem**2 + eta**2)

              end do
            end do

          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block D of the spinorbital manifold

  if(ispin == 4) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = - eGF(e) + eGF(m)
                num =     ERI(i,e,k,m)*ERI(j,m,l,e) -     ERI(i,e,k,m)*ERI(j,m,e,l) &
                    -     ERI(i,e,m,k)*ERI(j,m,l,e) +     ERI(i,e,m,k)*ERI(j,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) + 2d0*num*dem/(dem**2 + eta**2)

                dem = - eGF(e) + eGF(m)
                num =     ERI(j,e,k,m)*ERI(i,m,l,e) -     ERI(j,e,k,m)*ERI(i,m,e,l) &
                    -     ERI(j,e,m,k)*ERI(i,m,l,e) +     ERI(j,e,m,k)*ERI(i,m,e,l)

                KD_sta(ij,kl) = KD_sta(ij,kl) - 2d0*num*dem/(dem**2 + eta**2)

              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
