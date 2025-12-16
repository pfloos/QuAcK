subroutine compute_F_gen(O,V,N,Nsq,h,ERI_MO,rdm1,rdm2,F_gen)
      
! Compute generalized Fock-operator based on 1 and 2 rdm

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

  integer                       :: nind,m,p,q,r,s

  logical,parameter             :: debug = .false.

! Output variables

  double precision,intent(out)  :: F_gen(N,N)

! Compute Fgen

  F_gen(:,:) = 0d0
  
  do m=1,N
    do nind=1,N

      do q=1,N
        F_gen(m,nind) = F_gen(m,nind) + rdm1(m,q)*h(nind,q)
      end do

      do q=1,N
        do r=1,N
          do s=1,N
            F_gen(m,nind) = F_gen(m,nind) + rdm2(m,q,r,s)*ERI_MO(nind,r,q,s)
          end do
        end do
      end do

    end do
  end do

! Dump F_gen

  if(debug) then

    write(*,*) 'F_gen at the level:'
    call matout(N,N,F_gen)
    write(*,*)

  end if

end
