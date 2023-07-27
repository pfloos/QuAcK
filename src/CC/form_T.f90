subroutine form_T(nC,nO,nV,nR,delta_OOOVVV,ub,ubb,EcCCT)

! Compute (T) correction

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR

  double precision,intent(in)   :: delta_OOOVVV(nO,nO,nO,nV,nV,nV)
  double precision,intent(in)   :: ub(nO,nO,nO,nV,nV,nV)
  double precision,intent(in)   :: ubb(nO,nO,nO,nV,nV,nV)

! Local variables

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d

! Output variables

  double precision,intent(out)  :: EcCCT

  EcCCT = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do k=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR
            do c=1,nV-nR

              EcCCT = EcCCT                                &
                    + (ub(i,j,k,a,b,c) + ubb(i,j,k,a,b,c)) &
                    * ubb(i,j,k,a,b,c)/delta_OOOVVV(i,j,k,a,b,c) 

            end do
          end do
        end do
      end do
    end do
  end do

  EcCCT = - EcCCT/36d0

end subroutine 
