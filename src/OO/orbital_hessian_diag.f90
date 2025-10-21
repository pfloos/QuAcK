subroutine orbital_hessian_diag(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
      
! Compute the orbital hessian given 1- and 2- rdms

  implicit none
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  integer,intent(in)            :: Nsq
  double precision,intent(in)   :: h(N,N)
  double precision,intent(in)   :: ERI_MO(N,N,N,N)
  double precision,intent(in)   :: rdm1(N,N)
  double precision,intent(in)   :: rdm2(N,N,N,N)

! Local variables

  integer                       :: p,q,r,s,t,u,w
  integer                       :: pq

  logical,parameter             :: debug = .false.

  double precision,allocatable  :: tmp(:,:,:,:)

  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: hess(Nsq,1)

! Compute intermediate array

  allocate(tmp(N,N,N,N))

  tmp(:,:,:,:) = 0d0

  do p=1,N
    do q=1,N

      do r=1,N
        do s=1,N

          tmp(p,q,r,s) = - h(s,p)*rdm1(r,q) - h(q,r)*rdm1(p,s)
          
          do u=1,N

            tmp(p,q,r,s) = tmp(p,q,r,s) + 0.5d0*(                          &
                Kronecker_delta(q,r)*(h(u,p)*rdm1(u,s) + h(s,u)*rdm1(p,u)) &
              + Kronecker_delta(p,s)*(h(u,r)*rdm1(u,q) + h(q,u)*rdm1(r,u)) )

          end do

          do u=1,N
            do w=1,N

            tmp(p,q,r,s) = tmp(p,q,r,s) + ERI_MO(u,w,p,r)*rdm2(u,w,q,s) + ERI_MO(q,s,u,w)*rdm2(p,r,u,w)

            end do
          end do

          do t=1,N
            do u=1,N

            tmp(p,q,r,s) = tmp(p,q,r,s) - (                                   &
                ERI_MO(s,t,p,u)*rdm2(r,t,q,u) + ERI_MO(t,s,p,u)*rdm2(t,r,q,u) &
              + ERI_MO(q,u,r,t)*rdm2(p,u,s,t) + ERI_MO(q,u,t,r)*rdm2(p,u,t,s) )

            end do
          end do

          do t=1,N
            do u=1,N
              do w=1,N

                tmp(p,q,r,s) = tmp(p,q,r,s) + 0.5d0*(                                                    & 
                    Kronecker_delta(q,r)*(ERI_MO(u,w,p,t)*rdm2(u,w,s,t) + ERI_MO(s,t,u,w)*rdm2(p,t,u,w)) &
                  + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)) )

              end do
            end do
          end do

        end do
      end do

    end do
  end do

  ! Flatten Hessian matrix and add permutations

  pq = 0
  do p=1,N
    do q=1,N

      pq = pq + 1
      !hess(pq,1) = tmp(p,r,q,s) - tmp(r,p,q,s) - tmp(p,r,s,q) + tmp(r,p,s,q)  ! This element does not make any sense to me
      hess(pq,1) = tmp(p,q,p,q) - tmp(q,p,p,q) - tmp(p,q,q,p) + tmp(q,p,q,p)
    end do
  end do
end
