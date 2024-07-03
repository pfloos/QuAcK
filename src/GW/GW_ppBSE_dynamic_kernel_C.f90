subroutine GW_ppBSE_dynamic_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nS,nVV,lambda,eGW,Om,rho,OmBSE,KC_dyn,ZC_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d
  integer                       :: ab,cd

! Output variables

  double precision,intent(out)  :: KC_dyn(nVV,nVV)
  double precision,intent(out)  :: ZC_dyn(nVV,nVV)

! Initialization

  KC_dyn(:,:) = 0d0
  ZC_dyn(:,:) = 0d0

! Build dynamic A matrix

  if(ispin == 1) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a,nBas-nR
        ab = ab + 1
 
        cd = 0
        do c=nO+1,nBas-nR
          do d=c,nBas-nR
            cd = cd + 1
  
            do m=1,nS
              dem = OmBSE - eGW(c) - Om(m) - eGW(b)
              num = rho(a,c,m)*rho(b,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(c) - Om(m) - eGW(a)
              num = rho(b,c,m)*rho(a,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(d) - Om(m) - eGW(a)
              num = rho(a,c,m)*rho(b,d,m) 

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(d) - Om(m) - eGW(b)
              num = rho(b,c,m)*rho(a,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

           end do
           KC_dyn(ab,cd) = 2d0*KC_dyn(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
           ZC_dyn(ab,cd) = 2d0*ZC_dyn(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
          end do
        end do

      end do
    end do

  end if

 if(ispin == 2) then

    ab = 0
    do a=nO+1,nBas-nR
      do b=a+1,nBas-nR
        ab = ab + 1

        cd = 0
        do c=nO+1,nBas-nR
          do d=c+1,nBas-nR
            cd = cd + 1

            do m=1,nS

              dem = OmBSE - eGW(c) - Om(m) - eGW(b)
              num = rho(a,c,m)*rho(b,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(c) - Om(m) - eGW(a)
              num = rho(b,c,m)*rho(a,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) - num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(d) - Om(m) - eGW(a)
              num = rho(a,c,m)*rho(b,d,m) 

              KC_dyn(ab,cd) = KC_dyn(ab,cd) + num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) - num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

              dem = OmBSE - eGW(d) - Om(m) - eGW(b)
              num = rho(b,c,m)*rho(a,d,m)

              KC_dyn(ab,cd) = KC_dyn(ab,cd) - num*dem/(dem**2 + eta**2)
              ZC_dyn(ab,cd) = ZC_dyn(ab,cd) + num*(dem**2 - eta**2)/(dem**2 + eta**2)**2

            end do
            
            KC_dyn(ab,cd) = 2d0*KC_dyn(ab,cd)
            ZC_dyn(ab,cd) = 2d0*ZC_dyn(ab,cd)
          end do
        end do

      end do
    end do

  end if

end subroutine 
