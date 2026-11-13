subroutine GWiw_LSQ_setup(nfreqs,wcoord,Omin,Omax,nTest,Omtest,P)

! Pseudo-inverse used for the energy-specific quadrature along the imaginary axis (I. Duchemin and X. Blase).
!
! The imaginary-axis part of the contour-deformation self-energy needs, for every distance a = z - e_r,
! weights zk(a) such that   sum_k zk(a) [ -2 Om / (w_k^2 + Om^2) ]  =  1/(a +- Om)   for all Om in [Omin,Omax]
! (the poles of W). This routine builds P = pinv(B) with B(j,k) = -2 Om_j/(w_k^2+Om_j^2), so that zk = P t(a).

  implicit none

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nTest
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: Omin,Omax

! Output variables

  double precision,intent(out)  :: Omtest(nTest)
  double precision,intent(out)  :: P(nfreqs,nTest)

! Local variables

  integer                       :: j,k,m,info,lwork
  double precision,parameter    :: rc = 1d-12
  double precision,allocatable  :: B(:,:),U(:,:),VT(:,:),S(:),work(:)

  do j=1,nTest
    Omtest(j) = Omin*(Omax/Omin)**(dble(j-1)/dble(nTest-1))
  end do

  allocate(B(nTest,nfreqs),U(nTest,nfreqs),VT(nfreqs,nfreqs),S(nfreqs))

  do k=1,nfreqs
    do j=1,nTest
      B(j,k) = -2d0*Omtest(j)/(wcoord(k)**2+Omtest(j)**2)
    end do
  end do

  lwork = -1
  allocate(work(1))
  call dgesvd('S','S',nTest,nfreqs,B,nTest,S,U,nTest,VT,nfreqs,work,lwork,info)
  lwork = int(work(1))
  deallocate(work)
  allocate(work(lwork))
  call dgesvd('S','S',nTest,nfreqs,B,nTest,S,U,nTest,VT,nfreqs,work,lwork,info)
  if(info /= 0) then
    write(*,*) ' Error in the SVD of the quadrature matrix, info = ',info
    error stop
  end if

  P = 0d0
  do m=1,nfreqs
    if(S(m) > rc*S(1)) then
      do j=1,nTest
        do k=1,nfreqs
          P(k,j) = P(k,j) + VT(m,k)*U(j,m)/S(m)
        end do
      end do
    end if
  end do

  deallocate(B,U,VT,S,work)

end subroutine
