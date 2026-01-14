subroutine read_quadrature(nfreqs,ntimes,wcoord,wweight)

  integer                      :: kind_int

  integer,intent(in)           :: nfreqs
  integer,intent(out)          :: ntimes
  double precision,intent(out) :: wcoord(nfreqs)
  double precision,intent(out) :: wweight(nfreqs)

  ntimes = 0
  kind_int = 1
  call cgqf(nfreqs,kind_int,0d0,0d0,0d0,1d0,wcoord,wweight)
  wweight(:)=wweight(:)/((1d0-wcoord(:))**2d0)
  wcoord(:)=wcoord(:)/(1d0-wcoord(:))

end subroutine
