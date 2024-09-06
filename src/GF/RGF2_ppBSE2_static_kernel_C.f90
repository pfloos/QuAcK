subroutine RGF2_ppBSE2_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nVV,lambda,ERI,eGF,KC_sta)

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

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,c,d,e
  integer                       :: a0,aa,ab,cd

  double precision              :: eta2
  double precision, allocatable :: Om_tmp(:,:)
  
! Output variables

  double precision,intent(out)  :: KC_sta(nVV,nVV)

! Initialization

  KC_sta(:,:) = 0d0
  eta2 = eta * eta

  allocate(Om_tmp(nO,nV))

  ! Compute the energy differences and denominator once and store them in a temporary array
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(m,e,dem) SHARED(nC,nO,nBas,nR, eta2, eGF, Om_tmp)
  !$OMP DO
  do m=nC+1,nO
     do e=nO+1,nBas-nR
        dem = eGF(m) - eGF(e)
        Om_tmp(m,e) = dem / (dem*dem + eta2)
     enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
! Second-order correlation kernel for the block C of the singlet manifold

!  --- --- ---
!    OpenMP implementation
!  --- --- ---

  if(ispin == 1) then

    a0 = nBas - nR - nO
    !$OMP PARALLEL DEFAULT(NONE)                              &
    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, e, num) &
    !$OMP          SHARED(nO, nBas, nR, nC, a0, ERI, Om_tmp, KC_sta, lambda)
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
     
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) & 
                    -     ERI(a,m,e,c)*ERI(b,e,d,m) -     ERI(a,m,e,c)*ERI(b,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) + num * Om_tmp(m,e)
            
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) & 
                    -     ERI(b,m,e,c)*ERI(a,e,d,m) -     ERI(b,m,e,c)*ERI(a,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) + num * Om_tmp(m,e)
            
              end do
            end do

            KC_sta(ab,cd) = 2d0*lambda*KC_sta(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
            
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
!                 num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) & 
!                     -     ERI(a,m,e,c)*ERI(b,e,d,m) -     ERI(a,m,e,c)*ERI(b,e,m,d)
   
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
            
!                 dem = eGF(m) - eGF(e)
!                 num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) & 
!                     -     ERI(b,m,e,c)*ERI(a,e,d,m) -     ERI(b,m,e,c)*ERI(a,e,m,d)
   
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + num*dem/(dem**2 + eta**2)
            
!               end do
!            end do
   
!            KC_sta(ab,cd) = 2d0*lambda*KC_sta(ab,cd)/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
            
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
    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, m, e, num) &
    !$OMP          SHARED(nO, nBas, nR, nC, a0, ERI, Om_tmp, KC_sta)
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
     
                num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) &
                    -     ERI(a,m,e,c)*ERI(b,e,d,m) +     ERI(a,m,e,c)*ERI(b,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) + 2d0 * num * Om_tmp(m,e)
            
                num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) &
                    -     ERI(b,m,e,c)*ERI(a,e,d,m) +     ERI(b,m,e,c)*ERI(a,e,m,d)

                KC_sta(ab,cd) = KC_sta(ab,cd) - 2d0 * num * Om_tmp(m,e)
            
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
!   if(ispin == 2) then
   
!     ab = 0
!     do a=nO+1,nBas-nR
!       do b=a+1,nBas-nR
!         ab = ab + 1
   
!         cd = 0
!         do c=nO+1,nBas-nR
!           do d=c+1,nBas-nR
!             cd = cd + 1
   
!             do m=nC+1,nO
!               do e=nO+1,nBas-nR
     
!                 dem = eGF(m) - eGF(e)
!                 num = 2d0*ERI(a,m,c,e)*ERI(b,e,d,m) -     ERI(a,m,c,e)*ERI(b,e,m,d) &
!                     -     ERI(a,m,e,c)*ERI(b,e,d,m) +     ERI(a,m,e,c)*ERI(b,e,m,d)
   
!                 KC_sta(ab,cd) = KC_sta(ab,cd) + 2d0*num*dem/(dem**2 + eta**2)
            
!                 dem = eGF(m) - eGF(e)
!                 num = 2d0*ERI(b,m,c,e)*ERI(a,e,d,m) -     ERI(b,m,c,e)*ERI(a,e,m,d) &
!                     -     ERI(b,m,e,c)*ERI(a,e,d,m) +     ERI(b,m,e,c)*ERI(a,e,m,d)
   
!                 KC_sta(ab,cd) = KC_sta(ab,cd) - 2d0*num*dem/(dem**2 + eta**2)
            
!               end do
!             end do
   
!           end do
!         end do
   
!       end do
!     end do
   
  !   end if
  

  deallocate(Om_tmp)

end subroutine 
