subroutine GF2_ppBSE2_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nVV,lambda,ERI,eGF,KC_sta)

! Compute the resonant part of the static BSE2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d,e
  integer                       :: ab,cd

! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization

  KC_sta(:,:) = 0d0

! Second-order correlation kernel for the block C of the singlet manifold

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nBas-nR
          do d=c,nBas-nR
            cd = cd + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) & 
                    -     ERI(a,m,e,c)*ERI(b,e,d,m) + 2d0*ERI(a,m,e,c)*ERI(b,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
            
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) & 
                    -     ERI(b,m,e,c)*ERI(a,e,d,m) + 2d0*ERI(b,m,e,c)*ERI(a,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

! Second-order correlation kernel for the block C of the triplet manifold

  if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a+1,nBas-nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1
  
            do m=nC+1,nO
              do e=nO+1,nBas-nR
     
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) &
                    -     ERI(a,m,e,c)*ERI(b,e,d,m) +     ERI(a,m,e,c)*ERI(b,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) + 2d0*num*dem/(dem**2 + eta**2)
            
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) &
                    -     ERI(b,m,e,c)*ERI(a,e,d,m) +     ERI(b,m,e,c)*ERI(a,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) - 2d0*num*dem/(dem**2 + eta**2)
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
