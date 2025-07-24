subroutine RG0W0_rdm2(nOrb,nO,nS,lampl,rampl,lp,rp,lambda,t,rdm2)

! Compute 2-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: nOrb,nS,nO
double precision, intent(in)     :: lampl(nS,nOrb),rampl(nS,nOrb),rp(nOrb),lp(nOrb)

! Local
integer                          :: a,b,c,d,i,j,k,l
integer                          :: ia,jb,kc,ja,ib,ld

! Output
double precision,intent(out)     :: rdm2(nOrb,nOrb,nOrb,nOrb)

rdm2(:,:,:,:) = 0d0 

!! Occuppied occupied - virtual virtual
!do i=1,nO
!  do j=1,nO
!    do a=nO+1,nOrb
!      do b=nO+1,nOrb
!        ia = a - nO + i*nOrb 
!        jb = b - nO + j*nOrb 
!        rdm2(i,j,a,b) = 2*t(ia,jb)
!        do k=1,nO
!          do l=1,nO
!            do c=nO+1,nOrb
!              do d=nO+1,nOrb
!                kc= c - nO + k*nOrb 
!                ld= d - nO + l*nOrb
!                rdm2(i,j,a,b) = rdm2(i,j,a,b) + t(kc,ia)*t(jb,ld)*lambda(kc,ld)
!              end do
!            end do
!          end do
!        end do
!      end do
!    end do
!  end do
!end do
!
!! Virtual virtual occupied occupied
!do i=1,nO
!  do j=1,nO
!    do a=nO+1,nOrb
!      do b=nO+1,nOrb
!        ai = i + (a-nO)*nOrb 
!        bj = j + (b-nO)*nOrb 
!        !rdm2(a,b,i,j) = lambda(ai,bj)
!      end do
!    end do
!  end do
!end do
!
!! Virtual virtual occupied occupied
!do i=1,nO
!  do j=1,nO
!    do a=nO+1,nOrb
!      do b=nO+1,nOrb
!        do k=1,nO
!          do c=nO+1,nOrb
!            kc = c - nO + k*nOrb 
!            ja = a - nO + j*nOrb
!            ib = a - nO + j*nOrb
!            rdm2(i,a,b,j) = rdm2(i,a,b,j) + 0.5*lambda(kc,ja)*t(kc,ib) 
!          end do
!        end do
!      end do
!    end do
!  end do
!end do
end subroutine
