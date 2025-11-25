subroutine RG0W0_rdms_crpa(ispin,O,V,N,nS,X,Y,rdm1_rpa,rdm2_rpa)
! Compute RPA 1/2-Reduced-Density-Matrix based on RG0W0
  implicit none
  include 'parameters.h'

! Input
integer,intent(in)              :: ispin,O,V,nS,N
double precision,intent(in)     :: X(nS,nS),Y(nS,nS)

!Output
double precision,intent(out)    :: rdm1_rpa(N,N),rdm2_rpa(N,N,N,N)

!Local
double precision,allocatable    :: XXT(:,:),XYT(:,:),YYT(:,:)
integer                         :: i,a,jind,b,jb,ia

double precision,external     :: Kronecker_delta

allocate(XXT(nS,nS))
allocate(XYT(nS,nS))
allocate(YYT(nS,nS))

rdm2_rpa(:,:,:,:) = 0d0
rdm1_rpa(:,:)     = 0d0
YYT               = matmul(Y,transpose(Y)) 
XYT               = matmul(X,transpose(Y)) 
XXT               = matmul(X,transpose(X))

do i = 1, O
  do a = O+1, N
    do jind = 1, O
      do b = O+1, N
        jb = b - O + (jind - 1) * V 
        ia = a - O +    (i - 1) * V
        rdm2_rpa(b,i,jind,a) = rdm2_rpa(b,i,jind,a) &
                         + 2d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))
        rdm2_rpa(i,jind,a,b) = rdm2_rpa(i,jind,a,b) &
                         + 2d0*(XYT(jb,ia) + XYT(ia,jb))
        ! Contributions from fab*dij - fij*dab
        if(i==jind) then
          rdm1_rpa(a,b) = rdm1_rpa(a,b) &
                        + 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))
        endif
        if(a==b) then
          rdm1_rpa(jind,i) = rdm1_rpa(jind,i) &
                       - 0.5d0*(XXT(jb,ia) + YYT(jb,ia) - Kronecker_delta(jb,ia))
        endif
      enddo
    enddo
  enddo
enddo
deallocate(XXT,YYT,XYT)
end subroutine
