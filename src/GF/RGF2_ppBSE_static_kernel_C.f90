subroutine RGF2_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nVV,lambda,ERI,eGF,KC_sta)

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
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d,e
  integer                       :: a0,aa,ab,cd
  
! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization

  KC_sta(:,:) = 0d0
  
! Second-order correlation kernel for the block C of the singlet manifold

! --- --- ---
! OpenMP implementation
! --- --- ---

if(ispin == 1) then

  a0 = nBas - nR - nO
  !$OMP PARALLEL DEFAULT(NONE)                              &
  !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, e, num, dem) &
  !$OMP          SHARED(nO, nBas, nR, nC, a0, ERI, KC_sta, lambda, eGF, eta)
  !$OMP DO
  do a=nO+1,nBas-nR
    aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
    do b=a,nBas-nR
      ab = aa + b

      cd = 0
      do c=nO+1,nBas-nR
        do d=c,nBas-nR
          cd = cd + 1

          do m=nC+1,nO
            do e=nO+1,nBas-nR
   
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(a,m,c,e)*ERI(e,b,m,d) - ERI(a,m,c,e)*ERI(e,b,d,m)  &
                     -    ERI(a,m,e,c)*ERI(e,b,m,d) - ERI(a,m,e,c)*ERI(e,b,d,m)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                num = 2d0*ERI(a,e,c,m)*ERI(m,b,e,d) - ERI(a,e,c,m)*ERI(m,b,d,e)  &
                     -    ERI(a,e,m,c)*ERI(m,b,e,d) - ERI(a,e,m,c)*ERI(m,b,d,e)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
                                                                                 
                num = 2d0*ERI(b,m,c,e)*ERI(e,a,m,d) - ERI(b,m,c,e)*ERI(e,a,d,m)  &
                     -    ERI(b,m,e,c)*ERI(e,a,m,d) - ERI(b,m,e,c)*ERI(e,a,d,m)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                num = 2d0*ERI(b,e,c,m)*ERI(m,a,e,d) - ERI(b,e,c,m)*ERI(m,a,d,e)  &
                     -    ERI(b,e,m,c)*ERI(m,a,e,d) - ERI(b,e,m,c)*ERI(m,a,d,e)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
          
            end do
          end do

          KC_sta(ab,cd) = lambda*KC_sta(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
          
         end do
       end do

    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end if

!  --- --- ---
!  Naive implementation
!  --- --- ---

!   if(ispin == 1) then
   
!     ab = 0
!     do a=nO+1,nBas-nR
!       do b=a,nBas-nR
!         ab = ab + 1
   
!         cd = 0
!         do c=nO+1,nBas-nR
!           do d=c,nBas-nR
!             cd = cd + 1
   
!             do m=nC+1,nO
!               do e=nO+1,nBas-nR
   
!                 dem = eGF(m) - eGF(e)
!                 num = 2d0*ERI(a,m,c,e)*ERI(e,b,m,d) - ERI(a,m,c,e)*ERI(e,b,d,m)  &
!                      -    ERI(a,m,e,c)*ERI(e,b,m,d) - ERI(a,m,e,c)*ERI(e,b,d,m)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
   
!                 num = 2d0*ERI(a,e,c,m)*ERI(m,b,e,d) - ERI(a,e,c,m)*ERI(m,b,d,e)  &
!                      -    ERI(a,e,m,c)*ERI(m,b,e,d) - ERI(a,e,m,c)*ERI(m,b,d,e)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
                                                                                 
!                 num = 2d0*ERI(b,m,c,e)*ERI(e,a,m,d) - ERI(b,m,c,e)*ERI(e,a,d,m)  &
!                      -    ERI(b,m,e,c)*ERI(e,a,m,d) - ERI(b,m,e,c)*ERI(e,a,d,m)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
   
!                 num = 2d0*ERI(b,e,c,m)*ERI(m,a,e,d) - ERI(b,e,c,m)*ERI(m,a,d,e)  &
!                      -    ERI(b,e,m,c)*ERI(m,a,e,d) - ERI(b,e,m,c)*ERI(m,a,d,e)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
          
!               end do
!            end do
   
!            KC_sta(ab,cd) = lambda*KC_sta(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
          
!           end do
!         end do
   
!       end do
!     end do
   
!   end if

! Second-order correlation kernel for the block C of the triplet manifold

!  --- --- ---
!  OpenMP implementation
!  --- --- ---
  
if(ispin == 2) then

  a0 = nBas - nR - nO - 1
  !$OMP PARALLEL DEFAULT(NONE)                              &
  !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, e, num, dem) &
  !$OMP          SHARED(nO, nBas, nR, nC, a0, ERI, KC_sta, eGF, eta)
  !$OMP DO
  do a = nO+1, nBas-nR
    aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO - 1
    do b = a+1, nBas-nR
      ab = aa + b

      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1

          do m=nC+1,nO
            do e=nO+1,nBas-nR
   
                dem = eGF(m) - eGF(e)
                num = 2d0*ERI(a,m,c,e)*ERI(e,b,m,d) - ERI(a,m,c,e)*ERI(e,b,d,m)  &
                     -    ERI(a,m,e,c)*ERI(e,b,m,d) + ERI(a,m,e,c)*ERI(e,b,d,m)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)

                num = 2d0*ERI(a,e,c,m)*ERI(m,b,e,d) - ERI(a,e,c,m)*ERI(m,b,d,e)  &
                     -    ERI(a,e,m,c)*ERI(m,b,e,d) + ERI(a,e,m,c)*ERI(m,b,d,e)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
                                                                                 
                num = 2d0*ERI(b,m,c,e)*ERI(e,a,m,d) - ERI(b,m,c,e)*ERI(e,a,d,m)  &
                     -    ERI(b,m,e,c)*ERI(e,a,m,d) + ERI(b,m,e,c)*ERI(e,a,d,m)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)

                num = 2d0*ERI(b,e,c,m)*ERI(m,a,e,d) - ERI(b,e,c,m)*ERI(m,a,d,e)  &
                     -    ERI(b,e,m,c)*ERI(m,a,e,d) + ERI(b,e,m,c)*ERI(m,a,d,e)
                                                                                 
                KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2) 
          
            end do
          end do

        end do
      end do

    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end if

!  --- --- ---
!  Naive implementation
!  --- --- ---

! if(ispin == 2) then

!    ab = 0
!    do a=nO+1,nBas-nR
!      do b=a+1,nBas-nR
!        ab = ab + 1
  
!        cd = 0
!        do c=nO+1,nBas-nR
!          do d=c+1,nBas-nR
!            cd = cd + 1
  
!            do m=nC+1,nO
!              do e=nO+1,nBas-nR
  
!                 dem = eGF(m) - eGF(e)
!                 num = 2d0*ERI(a,m,c,e)*ERI(e,b,m,d) - ERI(a,m,c,e)*ERI(e,b,d,m)  &
!                      -    ERI(a,m,e,c)*ERI(e,b,m,d) + ERI(a,m,e,c)*ERI(e,b,d,m)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
  
!                 num = 2d0*ERI(a,e,c,m)*ERI(m,b,e,d) - ERI(a,e,c,m)*ERI(m,b,d,e)  &
!                      -    ERI(a,e,m,c)*ERI(m,b,e,d) + ERI(a,e,m,c)*ERI(m,b,d,e)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2) 
                                                                                 
!                 num = 2d0*ERI(b,m,c,e)*ERI(e,a,m,d) - ERI(b,m,c,e)*ERI(e,a,d,m)  &
!                      -    ERI(b,m,e,c)*ERI(e,a,m,d) + ERI(b,m,e,c)*ERI(e,a,d,m)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2)
  
!                 num = 2d0*ERI(b,e,c,m)*ERI(m,a,e,d) - ERI(b,e,c,m)*ERI(m,a,d,e)  &
!                      -    ERI(b,e,m,c)*ERI(m,a,e,d) + ERI(b,e,m,c)*ERI(m,a,d,e)
                                                                                 
!                 KC_sta(ab,cd) = KC_sta(ab,cd) - num*dem/(dem**2 + eta**2) 
         
!              end do
!            end do
  
!          end do
!        end do
  
!      end do
!    end do
  
!   end if
  
  if(ispin == 4) then
 
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
                num =       (ERI(a,m,c,e) - ERI(a,m,e,c)) * (ERI(e,b,m,d) - ERI(e,b,d,m))
                num = num + (ERI(a,e,c,m) - ERI(a,e,m,c)) * (ERI(m,b,e,d) - ERI(m,b,d,e))
                num = num - (ERI(b,m,c,e) - ERI(b,m,e,c)) * (ERI(e,a,m,d) - ERI(e,a,d,m))
                num = num - (ERI(b,e,c,m) - ERI(b,e,m,c)) * (ERI(m,a,e,d) - ERI(m,a,d,e))
                                                                              
                KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
          
              end do
            end do
 
          end do
        end do
 
      end do
    end do

  end if
  
! Second-order correlation kernel for the block C of the spinorbital manifold

  

! deallocate(Om_tmp)

end subroutine 
