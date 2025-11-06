subroutine GG0W0_rdm2_rpa(O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1_hf,rdm1_rpa,rdm2)

! Compute RPA 2-Reduced-Density-Matrix based in RG0W0
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)               :: N,nS,O,V
double precision, intent(in)     :: lampl(nS,N),rampl(nS,N),rp(N),lp(N)
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)
double precision, intent(in)     :: rdm1_hf(N,N),rdm1_rpa(N,N)

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
                kc = c - O + (k-1)*V 
                ld = d - O + (l-1)*V
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

! occupied virtual virtual occupied
do i=1,O
  do j=1,O
    do a=O+1,N
      do b=O+1,N
        do k=1,O
          do c=O+1,N
            kc = c - O + (k-1)*V
            ja = a - O + (j-1)*V
            ib = b - O + (i-1)*V
            rdm2(i,a,b,j) = rdm2(i,a,b,j) + 0.5d0*lambda(kc,ja)*t(kc,ib) 
          end do
        end do
      end do
    end do
  end do
end do

!! Contribution from 1rdm
!do p=1,N
!  do q=1,N
!    do r=1,N
!      do s=1,N
!        rdm2(p,q,r,s) = rdm2(p,q,r,s)                                           &
!                      + rdm1_rpa(p,r)*rdm1_hf(q,s) + rdm1_hf(p,r)*rdm1_rpa(q,s) & 
!                      - rdm1_rpa(p,s)*rdm1_hf(q,r) - rdm1_hf(p,s)*rdm1_rpa(q,r) 
!      enddo
!    enddo
!  enddo
!enddo
end subroutine
