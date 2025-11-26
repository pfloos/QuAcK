subroutine RG0W0_rdm1_rpax_triplet(O,V,N,nS,lambda,t,rdm1)

! Compute RPA 1-Reduced-Density-Matrix based on RG0W0
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)               :: N,nS,O,V
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)

! Local
integer                          :: a,b,i,j,c,k,d,l,kc,id,jd,lb,la,ia,ja,ib
integer                          :: nn

! Output
double precision,intent(out)     :: rdm1(N,N)

rdm1(:,:) = 0d0

! Occupied-occupied block
do i = 1,O
  do j = 1,O
    do a = 1,V
      ia = a + (i-1)*V
      ja = a + (j-1)*V
      rdm1(i,j) = rdm1(i,j) + t(ia,ja)*lambda(ia,ja)
    end do
  end do
end do

! Virtual-virtual block
do a = 1,V
  do b = 1,V
    do i = 1,O
      ia = a + (i-1)*V
      ib = b + (i-1)*V
      rdm1(O+a,O+b) = rdm1(O+a,O+b) - t(ia,ib)*lambda(ia,ib)
    end do
  end do
end do
end subroutine
