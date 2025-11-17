subroutine GG0W0_rdm1_mu(mu,O,V,N,nS,lampl,rampl,lp,rp,lambda,t,rdm1)

! Compute RPA 1-Reduced-Density-Matrix based in GG0W0
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)               :: mu,N,nS,O,V
double precision, intent(in)     :: lampl(N,N,N),rampl(N,N,N),rp(N),lp(N)
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS)

! Local
integer                          :: a,b,i,j,c,k,d,l,kc,id,jd,lb,la,ic,jc
integer                          :: e,m,f,g,nind,me,jg,ig,nf,oind,oa,ob
integer                          :: nn
double precision                 :: temp

! Output
double precision,intent(out)     :: rdm1(N,N)

!-----------------------!
! Correlation part      !
!-----------------------!

rdm1(:,:) = 0d0

! Occupied - Occupied
do i=1,O
  do j=1,O
    
    rdm1(i,j) = rdm1(i,j) - rp(i)*lp(j)
    
    do c=O+1,N
      do k=1,O
        rdm1(i,j) = rdm1(i,j) - rampl(c,k,i)*lampl(c,k,j) + rampl(c,i,k)*lampl(c,j,k)
      end do
    end do
    
    do c=O+1,N
      do d=O+1,N
        id = d - O + (i-1)*V
        jd = d - O + (j-1)*V
        rdm1(i,j) = rdm1(i,j) + lampl(i,d,c)*rampl(j,d,c)
      end do
    end do
    
    do e=O+1,N
      do m=1,O
        do f=O+1,N
          do nind=1,O
            do g=O+1,N
              me = e - O + (m - 1)*V
              nf = f - O + (m - 1)*V
              jg = g - O + (j - 1)*V
              ig = g - O + (i - 1)*V
              temp = 0d0
              do l=1,O
                temp = temp + lampl(f,nind,l)*rampl(e,m,l)
              end do
              do d=O+1,N
                temp = temp - rampl(nind,f,d)*lampl(m,e,d)
              end do
              rdm1(i,j) = rdm1(i,j) - lambda(me,jg)*t(ig,nf)*temp
            end do
          end do
        end do
      end do
    end do

    do e=O+1,N
      do m=1,O
        do f=O+1,N
          do nind=1,O
            do g=O+1,N
              me = e - O + (m - 1)*V
              nf = f - O + (nind - 1)*V
              jg = g - O + (j - 1)*V
              ig = g - O + (i - 1)*V
              temp = 0d0
              do l=1,O
                temp = temp + lampl(g,j,l)*rampl(e,m,l)
              end do
              do d=O+1,N
                temp = temp - rampl(j,g,d)*lampl(m,e,d)
              end do
              rdm1(i,j) = rdm1(i,j) - lambda(me,nf)*t(nf,ig)*temp
            end do
          end do
        end do
      end do
    end do
    
  end do
end do

! Virtual -- virtual
do a=O+1,N
  do b=O+1,N

    rdm1(a,b) = rdm1(a,b) - rp(a)*lp(b)

    do c=O+1,N
      do k=1,O
        rdm1(a,b) = rdm1(a,b) + rampl(k,c,b)*lampl(k,c,a) - rampl(k,b,c)*lampl(k,a,c)
      end do
    end do

    do k=1,O
      do l=1,O
        rdm1(a,b) = rdm1(a,b) + lampl(a,j,k)*rampl(b,l,k)
      end do
    end do

    do m=1,O
      do e=O+1,N
        do oind=1,O
          do nind=1,O
            do f=O+1,N
              me = e - O + (m - 1)*V
              oa = a - O + (oind - 1)*V
              ob = b - O + (oind - 1)*V
              nf = f - O + (nind - 1)*V
              temp = 0d0
              do l=1,O
                temp = temp + lampl(f,nind,l)*rampl(e,m,l)
              end do
              do d=O+1,N
                temp = temp - rampl(nind,f,d)*lampl(m,e,d)
              end do
              rdm1(a,b) = rdm1(a,b) - temp*lambda(me,oa)*t(ob,nf)
            end do
          end do
        end do
      end do
    end do

    do m=1,O
      do e=O+1,N
        do oind=1,O
          do nind=1,O
            do f=O+1,N
              me = e - O + (m - 1)*V
              oa = a - O + (oind - 1)*V
              ob = b - O + (oind - 1)*V
              nf = f - O + (nind - 1)*V
              temp = 0d0
              do l=1,O
                temp = temp + lampl(a,oind,l)*rampl(e,m,l)
              end do
              do d=O+1,N
                temp = temp - rampl(oind,a,d)*lampl(m,e,d)
              end do
              rdm1(a,b) = rdm1(a,b) - temp*lambda(me,nf)*t(nf,ob)
            end do
          end do
        end do
      end do
    end do
  
  end do
end do

do a=O+1,N
  do i=1,O
    rdm1(a,i) = rdm1(a,i) + rp(a)*lp(i)
    rdm1(i,a) = rdm1(i,a) + rp(i)*lp(a)
  end do
end do

!!! TODO: Contribution from xi
end subroutine
