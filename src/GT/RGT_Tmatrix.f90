subroutine RGT_Tmatrix(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,lambda,ERI,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,T)

! Compute the T-matrix tensor elements

  implicit none
  include 'parameters.h'

  ! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOOs, nOOt
  integer,intent(in)            :: nVVs, nVVt
  integer,intent(in)            :: isp_T
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Om1s(nVVs)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs)
  double precision,intent(in)   :: Om2s(nOOs)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs)
  double precision,intent(in)   :: Om1t(nVVt)
  double precision,intent(in)   :: rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: Om2t(nOOt)
  double precision,intent(in)   :: rho2t(nBas,nBas,nOOt)

  ! Local variables
  
  double precision,external     :: Kronecker_delta
  integer                       :: p,q,r,s
  integer                       :: c,d,k,l
  integer                       :: kl,cd
  
  ! Output variables

  double precision,intent(out)  :: T(nBas,nBas,nBas,nBas)

  ! Initialization

  T(:,:,:,:) = 0d0

  ! Elements aaaa

  if(isp_T == 1) then
     
     !$OMP PARALLEL &
     !$OMP SHARED(nC,nO,nBas,nR,T,ERI,nOOt,nVVt,rho1t,rho2t,Om1t,Om2t) &
     !$OMP PRIVATE(p,q,r,s,cd,kl) &
     !$OMP DEFAULT(NONE)
     !$OMP DO
     do s=nC+1,nBas-nR
        do r=nC+1,nBas-nR
           do q=nC+1,nBas-nR
              do p=nC+1,nBas-nR

                 T(p,q,r,s) = ERI(p,q,r,s) - ERI(p,q,s,r)

                 do cd=1,nVVt
                    T(p,q,r,s) = T(p,q,r,s) - rho1t(p,q,cd) * rho1t(r,s,cd) / Om1t(cd)
                 end do ! cd

                 do kl=1,nOOt
                    T(p,q,r,s) = T(p,q,r,s) + rho2t(p,q,kl) * rho2t(r,s,kl) / Om2t(kl)
                 enddo ! kl
           
              enddo ! p 
           enddo ! q
        enddo ! r
     enddo ! s
     !$OMP END DO
     !$OMP END PARALLEL

  endif 

  ! Elements abab 

  if(isp_T == 2) then 
     !$OMP PARALLEL &
     !$OMP SHARED(nC,nO,nBas,nR,T,ERI,nOOs,nOOt,nVVs,nVVt,rho1s,rho2s,Om1s,Om2s,rho1t,rho2t,Om1t,Om2t) &
     !$OMP PRIVATE(p,q,r,s,cd,kl) &
     !$OMP DEFAULT(NONE)
     !$OMP DO
     do s=nC+1,nBas-nR
        do r=nC+1,nBas-nR
           do q=nC+1,nBas-nR
              do p=nC+1,nBas-nR

                 T(p,q,r,s) = ERI(p,q,r,s)
                 
                 do cd=1,nVVs
                    T(p,q,r,s) = T(p,q,r,s) - 0.5d0 * rho1s(p,q,cd) * rho1s(r,s,cd) / Om1s(cd)
                 end do ! cd

                 do cd=1,nVVt
                    T(p,q,r,s) = T(p,q,r,s) - 0.5d0 * rho1t(p,q,cd) * rho1t(r,s,cd) / Om1t(cd)
                 end do ! cd
                 
                 do kl=1,nOOs
                    T(p,q,r,s) = T(p,q,r,s) + 0.5d0 * rho2s(p,q,kl) * rho2s(r,s,kl) / Om2s(kl)
                 enddo ! kl

                 do kl=1,nOOt
                    T(p,q,r,s) = T(p,q,r,s) + 0.5d0 * rho2t(p,q,kl) * rho2t(r,s,kl) / Om2t(kl)
                 enddo ! kl
           
              enddo ! p 
           enddo ! q
        enddo ! r
     enddo ! s
     !$OMP END DO
     !$OMP END PARALLEL

  endif

  ! Elements baab

  if(isp_T == 3) then 
     !$OMP PARALLEL &
     !$OMP SHARED(nC,nO,nBas,nR,T,ERI,nOOs,nOOt,nVVs,nVVt,rho1s,rho2s,Om1s,Om2s,rho1t,rho2t,Om1t,Om2t) &
     !$OMP PRIVATE(p,q,r,s,cd,kl) &
     !$OMP DEFAULT(NONE)
     !$OMP DO
     do s=nC+1,nBas-nR
        do r=nC+1,nBas-nR
           do q=nC+1,nBas-nR
              do p=nC+1,nBas-nR

                 T(p,q,r,s) = - ERI(p,q,s,r)
                 
                 do cd=1,nVVs
                    T(p,q,r,s) = T(p,q,r,s) + 0.5d0 * rho1s(p,q,cd) * rho1s(r,s,cd) / Om1s(cd)
                 end do ! cd

                 do cd=1,nVVt
                    T(p,q,r,s) = T(p,q,r,s) - 0.5d0 * rho1t(p,q,cd) * rho1t(r,s,cd) / Om1t(cd)
                 end do ! cd
                 
                 do kl=1,nOOs
                    T(p,q,r,s) = T(p,q,r,s) - 0.5d0 * rho2s(p,q,kl) * rho2s(r,s,kl) / Om2s(kl)
                 enddo ! kl

                 do kl=1,nOOt
                    T(p,q,r,s) = T(p,q,r,s) + 0.5d0 * rho2t(p,q,kl) * rho2t(r,s,kl) / Om2t(kl)
                 enddo ! kl
           
              enddo ! p 
           enddo ! q
        enddo ! r
     enddo ! s
     !$OMP END DO
     !$OMP END PARALLEL

  endif 
  
end subroutine
