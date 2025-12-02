subroutine RG0W0_transform_rdms_f_h(O,V,N,nS,rdm1,rdm2)
  
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)              :: O,V,N,nS
double precision,intent(in)     :: rdm1(N,N)

!Output
double precision,intent(inout)  :: rdm2(N,N,N,N)

!Local
integer                         :: p,q,r,s
double precision,allocatable    :: rdm1_hf(:,:)

allocate(rdm1_hf(N,N))

rdm1_hf(:,:) = 0d0
do p=1,O
  rdm1_hf(p,p) = 2d0
enddo

do p=1,N
  do q=1,N
    do r=1,N
      do s=1,N
        rdm2(p,q,r,s) = rdm2(p,q,r,s)               & 
                      + 2*rdm1(p,r) * rdm1_hf(q,s)  &
                      - rdm1(p,s) * rdm1_hf(q,r)
      end do
    end do
  end do
end do

deallocate(rdm1_hf)
end subroutine
