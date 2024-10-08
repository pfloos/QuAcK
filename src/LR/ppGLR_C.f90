subroutine ppGLR_C(nBas,nC,nO,nV,nR,nVV,lambda,e,ERI,Cpp)

! Compute the C matrix of the pp channel

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas) 
  
! Local variables

  double precision              :: eF
  double precision,external     :: Kronecker_delta

  integer                       :: a,b,c,d,ab,cd
  integer                       :: a0, aa
  double precision              :: e_ab, tmp_ab, delta_ac, tmp_cd

! Output variables

  double precision,intent(out)  :: Cpp(nVV,nVV)

! Define the chemical potential
 
! eF = e(nO) + e(nO+1)
  eF = 0d0

! Build C matrix for the singlet manifold

  !$OMP PARALLEL &
  !$OMP SHARED(Cpp,lambda,ERI,e,eF,nC,nO,nBas,nR) &
  !$OMP PRIVATE(c,d,a,b,ab,cd) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do c=nO+1,nBas-nR
    do d=c+1,nBas-nR
      cd = (c-(nO+1))*(nBas-nR-(nO+1)) - (c-1-(nO+1))*(c-(nO+1))/2 + d - c
        do a=nO+1,nBas-nR
          do b=a+1,nBas-nR
            ab = (a-(nO+1))*(nBas-nR-(nO+1)) - (a-1-(nO+1))*(a-(nO+1))/2 + b - a

            Cpp(ab,cd) = + (e(a) + e(b) - eF)*Kronecker_delta(a,c)*Kronecker_delta(b,d) & 
                 + lambda*(ERI(a,b,c,d) - ERI(a,b,d,c))

        end do
      end do
    end do
 end do
 !$OMP END DO
 !$OMP END PARALLEL

end subroutine 
