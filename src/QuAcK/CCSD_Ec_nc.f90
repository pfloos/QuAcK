subroutine CCSD_Ec_nc(nO,nV,t1,t2,Fov,OOVV,EcCCSD)

! Compute the CCSD correlatio energy in non-conanical form

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(in)   :: Fov(nO,nV)

  double precision,intent(in)   :: OOVV(nO,nO,nV,nV)

! Local variables

  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: EcCCSD

! Compute CCSD correlation energy

  EcCCSD = 0d0

! Singles contribution

  do i=1,nO
      do a=1,nO
     
      EcCCSD = EcCCSD + Fov(i,a)*t1(i,a)

    end do
  end do

! Doubles contribution

  do i=1,nO
    do j=1,nO
      do a=1,nO
        do b=1,nO
     
          EcCCSD = EcCCSD                              & 
                 + 0.5d0*OOVV(i,j,a,b)*t1(i,a)*t1(j,b) &
                 + 0.25d0*OOVV(i,j,a,b)*t2(i,j,a,b)

        end do
      end do
    end do
  end do

end subroutine CCSD_Ec_nc
