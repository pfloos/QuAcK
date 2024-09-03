subroutine RGF2_ppBSE2_dynamic_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOO,lambda,ERI,eGF,OmBSE,KD_dyn,ZD_dyn)

! Compute the resonant part of the dynamic BSE2 matrix

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
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,l,m
  integer                       :: e
  integer                       :: ij,kl

! Output variables

  double precision,intent(out)  :: KD_dyn(nOO,nOO)
  double precision,intent(out)  :: ZD_dyn(nOO,nOO)

! Initialization

  KD_dyn(:,:) = 0d0
  ZD_dyn(:,:) = 0d0

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
     
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(i,e,k,m)*ERI(j,m,l,e) -     ERI(i,e,k,m)*ERI(j,m,e,l) & 
                    -     ERI(i,e,m,k)*ERI(j,m,l,e) + 2d0*ERI(i,e,m,k)*ERI(j,m,e,l)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(j,e,k,m)*ERI(i,m,l,e) -     ERI(j,e,k,m)*ERI(i,m,e,l) & 
                    -     ERI(j,e,m,k)*ERI(i,m,l,e) + 2d0*ERI(j,e,m,k)*ERI(i,m,e,l)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(i,m,k,e)*ERI(j,e,l,m) -     ERI(i,m,k,e)*ERI(j,e,m,l) & 
                    -     ERI(i,m,e,k)*ERI(j,e,l,m) + 2d0*ERI(i,m,e,k)*ERI(j,e,m,l)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(j,m,k,e)*ERI(i,e,l,m) -     ERI(j,m,k,e)*ERI(i,e,m,l) & 
                    -     ERI(j,m,e,k)*ERI(i,e,l,m) + 2d0*ERI(j,m,e,k)*ERI(i,e,m,l)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

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
     
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(i,e,k,m)*ERI(j,m,l,e) - ERI(i,e,k,m)*ERI(j,m,e,l) - ERI(i,e,m,k)*ERI(j,m,l,e) 

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(k) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(j,e,k,m)*ERI(i,m,l,e) - ERI(j,e,k,m)*ERI(i,m,e,l) - ERI(j,e,m,k)*ERI(i,m,l,e)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(i)
                num = 2d0*ERI(i,m,k,e)*ERI(j,e,l,m) - ERI(i,m,k,e)*ERI(j,e,m,l) - ERI(i,m,e,k)*ERI(j,e,l,m)

                KD_dyn(ij,kl) = KD_dyn(ij,kl) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = - OmBSE + eGF(l) - eGF(e) + eGF(m) + eGF(j)
                num = 2d0*ERI(j,m,k,e)*ERI(i,e,l,m) - ERI(j,m,k,e)*ERI(i,e,m,l) - ERI(j,m,e,k)*ERI(i,e,l,m) 

                KD_dyn(ij,kl) = KD_dyn(ij,kl) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZD_dyn(ij,kl) = ZD_dyn(ij,kl) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
