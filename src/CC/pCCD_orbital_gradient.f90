subroutine pCCD_orbital_gradient(O,V,N,h,ERI_MO,rdm1,rdm2,grad)
      
! Compute the orbital gradient at the pCCD level

  implicit none
      
! Input variables

  integer,intent(in)            :: O
  integer,intent(in)            :: V
  integer,intent(in)            :: N
  double precision,intent(in)   :: h(N,N)
  double precision,intent(in)   :: ERI_MO(N,N,N,N)
  double precision,intent(in)   :: rdm1(N,N)
  double precision,intent(in)   :: rdm2(N,N,N,N)

! Local variables

  integer                       :: p,q,r,s,t
  integer                       :: pq

  logical,parameter             :: debug = .false.

! Output variables

  double precision,intent(out)  :: grad(N**2)

! Compute gradient

  grad(:) = 0d0

  pq = 0
  do p=1,N
    do q=1,N

      pq = pq + 1

      do r=1,N
        grad(pq) = grad(pq) + h(r,p)*rdm1(r,q)  - h(q,r)*rdm1(p,r)
      end do

      do r=1,N
        do s=1,N
          do t=1,N
            grad(pq) = grad(pq) + (ERI_MO(r,s,p,t)*rdm2(r,s,q,t) - ERI_MO(q,t,r,s)*rdm2(p,t,r,s))
          end do
        end do
      end do

    end do
  end do

  if(debug) then

    write(*,*) 'Orbital gradient at the pCCD level:'
    call matout(N,N,grad)
    write(*,*)

  end if

end
