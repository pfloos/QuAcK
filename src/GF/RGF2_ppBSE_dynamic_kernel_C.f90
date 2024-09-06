subroutine RGF2_ppBSE_dynamic_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nVV,lambda,ERI,eGF,OmBSE,KC_dyn,ZC_dyn)

! Compute the resonant part of the dynamic BSE@GF2 matrix

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
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d,e
  integer                       :: ab,cd

! Output variables

  double precision,intent(out)  :: KC_dyn(nVV,nVV)
  double precision,intent(out)  :: ZC_dyn(nVV,nVV)

! Initialization

  KC_dyn(:,:) = 0d0
  ZC_dyn(:,:) = 0d0

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
     
                dem = OmBSE - eGF(c) + eGF(m) - eGF(e) - eGF(b)
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) & 
                    -     ERI(a,m,e,c)*ERI(b,e,d,m) + 2d0*ERI(a,m,e,c)*ERI(b,e,m,d)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(c) + eGF(m) - eGF(e) - eGF(a)
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) & 
                    -     ERI(b,m,e,c)*ERI(a,e,d,m) + 2d0*ERI(b,m,e,c)*ERI(a,e,m,d)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(d) + eGF(m) - eGF(e) - eGF(a)
                num = 2d0*ERI(a,e,c,m)*ERI(b,m,d,e) -     ERI(a,e,c,m)*ERI(b,m,e,d) & 
                    -     ERI(a,e,m,c)*ERI(b,m,d,e) + 2d0*ERI(a,e,m,c)*ERI(b,m,e,d)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(d) + eGF(m) - eGF(e) - eGF(b)
                num = 2d0*ERI(b,e,c,m)*ERI(a,m,d,e) -     ERI(b,e,c,m)*ERI(a,m,e,d) & 
                    -     ERI(b,e,m,c)*ERI(a,m,d,e) + 2d0*ERI(b,e,c,m)*ERI(a,m,e,d)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
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
     
                dem = OmBSE - eGF(c) + eGF(m) - eGF(e) - eGF(b)
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) - ERI(a,m,c,e)*ERI(b,e,m,d) - ERI(a,m,e,c)*ERI(b,e,d,m) 

                KC_dyn(ab,cd) = KC_dyn(ab,cd) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(c) + eGF(m) - eGF(e) - eGF(a)
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) - ERI(b,m,c,e)*ERI(a,e,m,d) - ERI(b,m,e,c)*ERI(a,e,d,m)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(d) + eGF(m) - eGF(e) - eGF(a)
                num = 2d0*ERI(a,e,c,m)*ERI(b,m,d,e) - ERI(a,e,c,m)*ERI(b,m,e,d) - ERI(a,e,m,c)*ERI(b,m,d,e)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) + 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
                dem = OmBSE - eGF(d) + eGF(m) - eGF(e) - eGF(b)
                num = 2d0*ERI(b,e,c,m)*ERI(a,m,d,e) - ERI(b,e,c,m)*ERI(a,m,e,d) - ERI(b,e,m,c)*ERI(a,m,d,e)

                KC_dyn(ab,cd) = KC_dyn(ab,cd) - 0.5d0*num*dem/(dem**2 + eta**2)
                ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + 0.5d0*num*(dem**2 - eta**2)/(dem**2 + eta**2)**2
            
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
