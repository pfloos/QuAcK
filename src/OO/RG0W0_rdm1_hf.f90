subroutine RG0W0_rdm1_hf(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1)

! Compute HF 1-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: N,nS,O,V
double precision, intent(in)     :: lampl(nS,N),rampl(nS,N),rp(N),lp(N)
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)

! Local
integer                          :: a,b,i,j,c,k,d,l,kc,id,jd,lb,la
integer                          :: nn

! Output
double precision,intent(out)     :: rdm1(N,N)

rdm1(:,:) = 0d0

! Occupied
do i=1,O
  rdm1(i,i) = rdm1(i,i) + 2d0 ! HF contribution
end do
end subroutine
