subroutine RG0W0_rdm2_rpax_triplet(O,V,N,nS,lambda,t,rdm2)

! Compute RPA 2-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: N,nS,O,V
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)

! Local
integer                          :: a,b,c,d,i,j,k,l,p,q,r,s
integer                          :: ia,jb,kc,ja,ib,ld
integer                          :: nn

! Output
double precision,intent(out)     :: rdm2(N,N,N,N)

rdm2(:,:,:,:) = 0d0

! Occuppied occupied - virtual virtual
do i=1,O
  do j=1,O
    do a=O+1,N
      do b=O+1,N
        ia = a - O + (i-1)*V 
        jb = b - O + (j-1)*V
        rdm2(i,j,a,b) = rdm2(i,j,a,b) + t(ia,jb)
        do k=1,O
          do l=1,O
            do c=O+1,N
              do d=O+1,N
                kc= c - O + (k-1)*V 
                ld= d - O + (l-1)*V
                rdm2(i,j,a,b) = rdm2(i,j,a,b) + t(kc,ia)*t(jb,ld)*lambda(kc,ld)
              end do
            end do
          end do
        end do
      end do
    end do
  end do
end do

! Virtual virtual occupied occupied
do i=1,O
  do j=1,O
    do a=O+1,N
      do b=O+1,N
        ia = a - O + (i-1)*V 
        jb = b - O + (j-1)*V
        rdm2(a,b,i,j) = rdm2(a,b,i,j) + lambda(ia,jb)
      end do
    end do
  end do
end do
end subroutine
