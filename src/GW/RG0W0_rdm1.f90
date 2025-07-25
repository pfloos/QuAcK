subroutine RG0W0_rdm1(nOrb,nO,nS,lampl,rampl,lp,rp,lambda,t,rdm1)

! Compute 1-Reduced-Density-Matrix based in RG0W0

! Input
integer,intent(in)               :: nOrb,nS,nO
double precision, intent(in)     :: lampl(nS,nOrb),rampl(nS,nOrb),rp(nOrb),lp(nOrb)
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)

! Local
integer                          :: a,b,i,j,c,k,d,l,kc,id,jd,lb,la

! Output
double precision,intent(out)     :: rdm1(nOrb,nOrb)

rdm1(:,:) = 0d0

!! Occupied
!do i=1,nO
!  do j=1,nO
!    do c=nO+1,nOrb
!      do k=1,nO
!        do d=nO+1,nOrb
!          kc = c - nO + k*nOrb
!          id = d - nO + i*nOrb
!          jd = d - nO + j*nOrb
!          rdm1(i,j) = rdm1(i,j) - 0.5*(t(kc,id)*lambda(kc,jd))
!        end do
!      end do
!    end do
!  end do
!end do
!
!! Virtual
!do a=nO+1,nOrb
!  do b=nO+1,nOrb
!    do c=nO+1,nOrb
!      do k=1,nO
!        do l=1,nOrb
!          kc = c - nO + k*nOrb
!          lb = b - nO + l*nOrb
!          la = a - nO + l*nOrb
!          rdm1(a,b) = rdm1(a,b) + 0.5*(t(kc,lb)*lambda(kc,la))
!        end do
!      end do
!    end do
!  end do
!end do
end subroutine
