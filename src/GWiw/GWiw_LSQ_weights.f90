subroutine GWiw_LSQ_weights(nfreqs,nTest,P,Omtest,a,zk)

! Energy-specific imaginary-axis quadrature weights zk(a) for the complex distance a = z - e_r
! (see GWiw_LSQ_setup). The target is 1/(a+Om) for Re(a) >= 0 and 1/(a-Om) for Re(a) < 0.

  implicit none

! Input variables

  integer,intent(in)            :: nfreqs,nTest
  double precision,intent(in)   :: P(nfreqs,nTest)
  double precision,intent(in)   :: Omtest(nTest)
  complex(kind=8),intent(in)    :: a

! Output variables

  complex(kind=8),intent(out)   :: zk(nfreqs)

! Local variables

  integer                       :: j,k
  double precision              :: sg
  complex(kind=8)               :: t

  sg = 1d0
  if(real(a) < 0d0) sg = -1d0

  zk = (0d0,0d0)
  do j=1,nTest
    t = 1d0/(a+sg*Omtest(j))
    do k=1,nfreqs
      zk(k) = zk(k) + P(k,j)*t
    end do
  end do

end subroutine
