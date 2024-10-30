subroutine GGT_Tmatrix(nOrb,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,Om1,rho1,Om2,rho2,T)

! Compute the T-matrix tensor elements

  implicit none
  include 'parameters.h'

  ! Input variables

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eGF(nOrb)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nOrb,nOrb,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nOrb,nOrb,nOO)

  ! Local variables

  integer                       :: p,q,r,s
  integer                       :: c,d,k,l
  integer                       :: kl,cd
  
  ! Output variables

  double precision,intent(out)  :: T(nOrb,nOrb,nOrb,nOrb)

  ! Initialization
  T(:,:,:,:) = 0d0
  
! Start by building the tensor elements of T
! This is probabbly not a good idea because this tensor is really large
  !$OMP PARALLEL &
  !$OMP SHARED(nC,nO,nOrb,nR,T,ERI,rho1,rho2,Om1,Om2) &
  !$OMP PRIVATE(p,q,r,s,c,d,cd,k,l,kl) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do s=nC+1,nOrb-nR
    do r=nC+1,nOrb-nR
      do q=nC+1,nOrb-nR
        do p=nC+1,nOrb-nR
          T(p,q,r,s) = ERI(p,q,r,s) - ERI(p,q,s,r)
           
          cd = 0
          do c = nO+1, nOrb-nR
            do d = c+1, nOrb-nR
              cd = cd + 1
              T(p,q,r,s) = T(p,q,r,s) - rho1(p,q,cd) * rho1(r,s,cd) / Om1(cd)
            end do ! d
          end do ! c

          kl = 0
          do k = nC+1, nO
            do l = k+1, nO
              kl = kl + 1
              T(p,q,r,s) = T(p,q,r,s) + rho2(p,q,kl) * rho2(r,s,kl) / Om2(kl)
            enddo
          enddo
           
        enddo
      enddo
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  
end subroutine
