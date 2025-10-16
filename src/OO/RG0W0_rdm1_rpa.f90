subroutine RG0W0_rdm1_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1)

! Compute RPA 1-Reduced-Density-Matrix based on RG0W0

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
  do j=1,O
    do c=O+1,N
      do k=1,O
        do d=O+1,N
          kc = c - O + (k-1)*V
          id = d - O + (i-1)*V
          jd = d - O + (j-1)*V
          rdm1(i,j) = rdm1(i,j) - 0.5d0*(t(kc,id)*lambda(kc,jd))
        end do
      end do
    end do
  end do
end do

! Virtual
do a=O+1,N
  do b=O+1,N
    do c=O+1,N
      do k=1,O
        do l=1,O
          kc = c - O + (k-1)*V
          lb = b - O + (l-1)*V
          la = a - O + (l-1)*V
          rdm1(a,b) = rdm1(a,b) + 0.5d0*(t(kc,lb)*lambda(kc,la))
        end do
      end do
    end do
  end do
end do
end subroutine
