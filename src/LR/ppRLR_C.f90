subroutine ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVV,lambda,e,ERI,Cpp)

! Compute the C matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nOrb),ERI(nOrb,nOrb,nOrb,nOrb) 
  
! Local variables

  double precision              :: eF
  double precision,external     :: Kronecker_delta

  logical                       :: use_chem_pot
  integer                       :: a,b,c,d,ab,cd
  integer                       :: a0, aa
  double precision              :: e_ab, tmp_ab, delta_ac, tmp_cd

! Output variables

  double precision,intent(out)  :: Cpp(nVV,nVV)

! Define the chemical potential
 
  use_chem_pot=.false.
  eF = 0d0
  inquire(file='use_chem_pot_in_ppRPA', exist=use_chem_pot)
  if(use_chem_pot) then
    write(*,*) 'File use_chem_pot_in_ppRPA encountered, setting eF = eHOMO + eLUMO'
    eF = e(nO) + e(nO+1)
  endif

! Build C matrix for the singlet manifold

  if(ispin == 1) then

    a0 = nOrb - nR - nO

    !$OMP PARALLEL DEFAULT(NONE)                                                   &
    !$OMP          PRIVATE(a, b, aa, ab, c, d, cd, e_ab, tmp_ab, delta_ac, tmp_cd) &
    !$OMP          SHARED(nO, nOrb, nR, a0, eF, lambda, e, ERI, Cpp)
    !$OMP DO
    do a = nO+1, nOrb-nR
      aa = a0 * (a - nO - 1) - (a - nO - 1) * (a - nO) / 2 - nO
      do b = a, nOrb-nR
        ab = aa + b

        e_ab = e(a) + e(b) - eF

        tmp_ab = lambda
        if(a .eq. b) then
          tmp_ab = 0.7071067811865475d0 * lambda
        endif

        cd = 0
        do c = nO+1, nOrb-nR

          delta_ac = 0.d0
          if(a .eq. c) then
            delta_ac = 1.d0
          endif

         do d = c, nOrb-nR
            cd = cd + 1

            tmp_cd = tmp_ab
            if(c .eq. d) then
              tmp_cd = 0.7071067811865475d0 * tmp_ab
            endif

            Cpp(ab,cd) = 0.d0
            if(b .eq. d) then
              Cpp(ab,cd) = e_ab * delta_ac
            endif
 
            Cpp(ab,cd) = Cpp(ab,cd) + tmp_cd * (ERI(a,b,c,d) + ERI(a,b,d,c))
          enddo
        enddo
      enddo
    enddo
    !$OMP END DO
    !$OMP END PARALLEL


!    ab = 0
!    do a=nO+1,nOrb-nR
!     do b=a,nOrb-nR
!        ab = ab + 1
!        cd = 0
!        do c=nO+1,nOrb-nR
!         do d=c,nOrb-nR
!            cd = cd + 1
! 
!            Cpp(ab,cd) = + (e(a) + e(b) - eF)*Kronecker_delta(a,c)*Kronecker_delta(b,d) &
!                         + lambda*(ERI(a,b,c,d) + ERI(a,b,d,c))/sqrt((1d0 + Kronecker_delta(a,b))*(1d0 + Kronecker_delta(c,d)))
! 
!          end do
!        end do
!      end do
!    end do

  end if

! Build C matrix for the triplet or alpha-alpha manifold

  if(ispin == 2 .or. ispin == 4) then
    !$OMP PARALLEL &
    !$OMP SHARED(Cpp,lambda,ERI,e,eF,nC,nO,nOrb,nR) &
    !$OMP PRIVATE(c,d,a,b,ab,cd) &
    !$OMP DEFAULT(NONE)
    !$OMP DO
    do c=nO+1,nOrb-nR
      do d=c+1,nOrb-nR
        cd = (c-(nO+1))*(nOrb-nR-(nO+1)) - (c-1-(nO+1))*(c-(nO+1))/2 + d - c
          do a=nO+1,nOrb-nR
            do b=a+1,nOrb-nR
              ab = (a-(nO+1))*(nOrb-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a
 
              Cpp(ab,cd) = + (e(a) + e(b) - eF)*Kronecker_delta(a,c)*Kronecker_delta(b,d) & 
                   + lambda*(ERI(a,b,c,d) - ERI(a,b,d,c))
 
          end do
        end do
      end do
   end do
   !$OMP END DO
   !$OMP END PARALLEL

 end if

! Build the alpha-beta block of the C matrix

  if(ispin == 3) then

    ab = 0
    do a=nO+1,nOrb-nR
     do b=nO+1,nOrb-nR
        ab = ab + 1
        cd = 0
        do c=nO+1,nOrb-nR
         do d=nO+1,nOrb-nR
            cd = cd + 1
 
            Cpp(ab,cd) = + (e(a) + e(b) - eF)*Kronecker_delta(a,c)*Kronecker_delta(b,d) & 
                         + lambda*ERI(a,b,c,d)
 
          end do
        end do
      end do
    end do

  end if

end subroutine 
