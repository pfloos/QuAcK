subroutine RG0W0_rdm2_hf(O,V,N,nS,rdm2)

! Compute HF 2-Reduced-Density-Matrix based in RG0W0

implicit none
include 'parameters.h'

! Input
integer,intent(in)               :: N,nS,O,V

! Local
integer                          :: a,b,c,d,i,j,k,l,p,q,r,s
integer                          :: ia,jb,kc,ja,ib,ld
integer                          :: nn

! Output
double precision,intent(out)     :: rdm2(N,N,N,N)

rdm2(:,:,:,:) = 0d0

! HF contribution
do p=1,O
  do q=1,O
    do r=1,O
      do s=1,O
        rdm2(p,q,r,s) = 2d0*DBLE(2*MERGE(1.0d0, 0.0d0, p==r) * MERGE(1.0d0,0.0d0,q==s) - &
                     MERGE(1.0d0,0.0d0,p==s) * MERGE(1.0d0,0.0d0,q==r))
        end do
    end do
  end do
end do
end subroutine
