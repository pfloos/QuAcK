subroutine GG0W0_rdm2_mu(O,V,N,nS,lampl,rampl,lp,rp,xi,lambda,t,rdm2)
   

! Compute RPA 2-Reduced-Density-Matrix based in RG0W0
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)               :: N,nS,O,V
double precision, intent(in)     :: lambda(nS,nS),t(nS,nS),xi(nS,nS)
double precision, intent(in)     :: lampl(N,N,N),rampl(N,N,N),rp(N),lp(N)

! Local
integer                          :: a,b,c,d,i,j,k,l,p,q,r,s,e,f,nind,m
integer                          :: ia,jb,kc,ja,ib,ld,nf,me
integer                          :: nn
double precision                 :: tmp

! Output
double precision,intent(out)     :: rdm2(N,N,N,N)

rdm2(:,:,:,:) = 0d0

do i=1,O
  do j=1,O
    do a=O+1,N
      do b=O+1,N

        ia = a - O + (i-1)*V 
        ib = b - O + (i-1)*V 
        jb = b - O + (j-1)*V
        ja = a - O + (j-1)*V
        
        ! Virtual-virtual Occupied-occupied
        rdm2(a,b,i,j) =  - rampl(i,a,b)*lp(j) & 
                         - lampl(a,i,j)*rp(b) 
        do c=O+1,N
          do k=1,O
            kc = c - O + (k-1)*V
            tmp = 0d0
            do l=1,O
              tmp = tmp + lampl(b,j,l)*rampl(c,k,l)
            end do
            do d=O+1,N
              tmp = tmp - rampl(j,b,d)*lampl(k,c,d) 
            end do
            rdm2(a,b,i,j) = rdm2(a,b,i,j) + lambda(kc,ia)*tmp
          end do
        end do

        ! Occupied-occupied Virtual -virtual
        rdm2(i,j,a,b) =  - rampl(a,i,j)*lp(b) &
                         - lampl(i,a,b)*rp(j)
        do c=O+1,N
          do k=1,O
            kc = c - O + (k-1)*V
            do d=O+1,N
              do l=1,O
                ld = d - O + (l-1)*V
                
                rdm2(i,j,a,b) =  rdm2(i,j,a,b)                              &
                                - t(jb,ld)*lambda(kc,ld)*rampl(c,k,i)*lp(a) &
                                - t(jb,ld)*lambda(kc,ld)*lampl(k,c,a)*rp(i)
              end do
            end do
            tmp = 0d0
            do l=1,O
              tmp = tmp + rampl(b,j,l)*lampl(c,k,l)
            end do
            do d=O+1,N
              tmp = tmp - lampl(j,b,d)*rampl(k,c,d) 
            end do
            rdm2(i,j,a,b) = rdm2(i,j,a,b) + t(kc,ia)*tmp  
          end do
        end do
        
        do c =O+1,N
          do e=O+1,N
            do f=O+1,N
              do k=1,O
                kc = c - O + (k-1)*V
                do m=1,O
                  me = e - O + (m-1)*V
                  do nind=1,O
                    
                    nf = f - O + (nind-1)*V
                    tmp = 0d0
                    do l=1,O
                      tmp = tmp + rampl(f,nind,l)*lampl(c,k,l)
                    end do
                    do d=O+1,N
                      tmp = tmp - lampl(nind,f,d)*rampl(k,c,d) 
                    end do
                    rdm2(i,j,a,b) = rdm2(i,j,a,b) + tmp*t(kc,ia)*t(jb,me)*lambda(me,nf) 

                  end do
                end do
              end do
            end do
          end do
        end do

        ! Occ-virt virt-occ
        do l=1,O
          rdm2(i,a,b,j) = rdm2(i,a,b,j) + 0.5d0*lampl(a,j,l)*rampl(b,i,l) 
        end do
        do d=O+1,N
          rdm2(i,a,b,j) = rdm2(i,a,b,j) - 0.5d0*rampl(j,a,d)*lampl(i,b,d)
        end do
        do d=O+1,N
          do l=1,O
            ld = d - O + (l-1)*V
            rdm2(i,a,b,j) = rdm2(i,a,b,j) - 0.5d0*lambda(ja,ld)*rampl(d,l,i)*lp(b)
            rdm2(i,a,b,j) = rdm2(i,a,b,j) - 0.5d0*lambda(ja,ld)*lampl(l,d,b)*rp(i)
            rdm2(i,a,b,j) = rdm2(i,a,b,j) - 0.5d0*t(ib,ld)*rampl(l,d,a)*lp(j)
            rdm2(i,a,b,j) = rdm2(i,a,b,j) - 0.5d0*t(ib,ld)*lampl(d,l,j)*rp(a)
          end do
        end do
        do c=O+1,N
          do e=O+1,N
            do k=1,O
              kc = c - O + (k-1)*V
              do m=1,O
                tmp = 0d0
                do l=1,O
                  tmp = tmp + lampl(a,j,l)*rampl(e,m,l)
                end do
                do d=O+1,N
                  tmp = tmp - rampl(j,a,d)*lampl(m,e,d) 
                end do
                rdm2(i,a,b,j) = rdm2(i,a,b,j) + 0.5d0*t(ib,kc)*lambda(kc,me)*tmp
                tmp = 0d0 
                do l=1,O
                  tmp = tmp + lampl(c,k,l)*rampl(e,m,l) 
                end do
                do d=O+1,N
                  tmp = tmp - rampl(k,c,d)*lampl(m,e,d) 
                end do
                rdm2(i,a,b,j) = rdm2(i,a,b,j) + 0.5d0*t(ib,kc)*lambda(ja,me)*tmp
              end do
            end do
          end do
        end do
      
      end do
    end do
  end do
end do

! P(ia,jb) on OOVV and VVOO Block
do i=1,O
  do j=1,O
    do a=O+1,N
      do b=O+1,N
        rdm2(i,j,a,b) = rdm2(i,j,a,b) + rdm2(j,i,b,a)
        rdm2(a,b,i,j) = rdm2(a,b,i,j) + rdm2(b,a,j,i)
      end do
    end do
  end do
end do


do i=1,O
  do j=1,O
    do k=1,O
      do a=O+1,N

        ja = a - O + (j-1)*V

        ! Occ Occ Occ Virt
        rdm2(i,j,k,a) = rdm2(i,j,k,a) - rampl(a,j,i)*lp(k)
        do e=O+1,N
          do m=1,O
            me = e - O + (m-1)*V
            rdm2(i,j,k,a) = rdm2(i,j,k,a) - t(me,ja)*rp(i)*lampl(e,m,k)
            do d=O+1,N
              do l=1,O
                ld = d - O + (l-1)*V
                rdm2(i,j,k,a) = rdm2(i,j,k,a) - t(me,ja)*lambda(me,ld)*rampl(d,l,i)*lp(k)
              end do
            end do
          end do
        end do
        
        ! Occ Virt Occ Occ
        rdm2(k,a,i,j) = rdm2(k,a,i,j) - lampl(a,j,i)*rp(k)
        do e=O+1,N
          do m=1,O
            me = e - O + (m-1)*V
            rdm2(k,a,i,j) = rdm2(k,a,i,j) - lambda(me,ja)*lp(i)*rampl(e,m,k)
          end do
        end do

      end do
    end do
  end do
end do 

do a=O+1,N
  do b=O+1,N
    do c=O+1,N
      do i=1,O
        ib = b - O + (i-1)*V
        ! Virt occ virt virt
        rdm2(c,i,a,b) = rdm2(c,i,a,b) - lampl(i,b,a)*rp(c)
        do e=O+1,N
          do m=1,O
            me = e - O + (m-1)*V
            rdm2(c,i,a,b) = rdm2(c,i,a,b) - t(me,ib)*lp(a)*rampl(m,e,c)
            do d=O+1,N
              do l=1,O
                ld = d - O + (l-1)*V
                rdm2(c,i,a,b) = rdm2(c,i,a,b) - t(me,ib)*lambda(me,ld)*lampl(l,d,a)*rp(c)
              end do
            end do
          end do
        end do

        ! Virt virt virt occ
        rdm2(a,b,c,i) = rdm2(a,b,c,i) ! continue here
        
      end do
    end do
  end do
end do
end subroutine
