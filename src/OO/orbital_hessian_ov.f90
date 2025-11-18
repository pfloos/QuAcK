subroutine orbital_hessian_ov(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,hess)
      
! Computes the ovov part of the orbital hessian given 1- and 2- rdms the remaining entries are set to 0

  implicit none
  include 'parameters.h'
      
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

  integer                       :: p,q,r,s,t,u,w,i,j,a,b
  integer                       :: pq,rs

  logical,parameter             :: debug = .false.

  double precision,allocatable  :: tmp(:,:,:,:)
  double precision,allocatable  :: hess2(:,:)
  double precision,allocatable  :: e(:)

  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: hess(Nsq,Nsq)

! Compute intermediate array

  allocate(tmp(N,N,N,N))

  tmp(:,:,:,:) = 0d0

  do i=1,O
    do a=O+1,N

      do j=1,O
        do b=O+1,N
          ! first permutation
          p = i
          q = a 
          r = j
          s = b
          
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
                  + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)))

              end do
            end do
          end do

          ! second permutation
          p = a
          q = i 
          r = j
          s = b
          
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
                  + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)))

              end do
            end do
          end do

          ! third permutation
          p = i
          q = a 
          r = b
          s = j
          
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
                  + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)))

              end do
            end do
          end do

          ! fourth permuation
          p = a
          q = i 
          r = b
          s = j
          
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
                  + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)))

              end do
            end do
          end do

        end do
      end do

    end do
  end do

  ! Flatten Hessian matrix and add permutations

  do p=1,O
    do q=O+1,N

      pq = q + (p-1)*N
   
      rs = 0
      do r=1,O
        do s=O+1,N

          rs = s + (r-1)*N

!         hess(pq,rs) = tmp(p,r,q,s) - tmp(r,p,q,s) - tmp(p,r,s,q) + tmp(r,p,s,q)
         hess(pq,rs) = tmp(p,q,r,s) - tmp(q,p,r,s) - tmp(p,q,s,r) + tmp(q,p,s,r)

        end do
      end do

    end do
  end do
deallocate(tmp)

! Dump Hessian

  if(debug) then

    write(*,*) 'Orbital Hessian at the level:'
    call matout(Nsq,Nsq,hess)
    write(*,*)

  end if

end
