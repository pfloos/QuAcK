subroutine form_r1_nc(nO,nV,t1,t2,delta_ov,Fov,cFoo,cFov,cFvv,OOVO,OVOV,OVVV,r1)

! Form residues for t1 in non-canonical CCSD

  implicit none

! Input variables

  integer,intent(in)            :: nO,nV

  double precision,intent(in)   :: t1(nO,nV)
  double precision,intent(in)   :: t2(nO,nO,nV,nV)

  double precision,intent(in)   :: delta_ov(nO,nV)
  double precision,intent(in)   :: Fov(nO,nV)

  double precision,intent(in)   :: cFoo(nO,nO)
  double precision,intent(in)   :: cFov(nO,nV)
  double precision,intent(in)   :: cFvv(nV,nV)

  double precision,intent(in)   :: OOVO(nO,nO,nV,nO)
  double precision,intent(in)   :: OVOV(nO,nV,nO,nV)
  double precision,intent(in)   :: OVVV(nO,nV,nV,nV)

! Local variables

  integer                       :: i,j,m,n
  integer                       :: a,b,e,f

! Output variables

  double precision,intent(out)  :: r1(nO,nV)

  r1(:,:) = Fov(:,:)

  do i=1,nO
    do a=1,nV

      do e=1,nV
        r1(i,a) = r1(i,a) + t1(i,e)*cFvv(a,e)
      end do

      do m=1,nO
        r1(i,a) = r1(i,a) - t1(m,a)*cFoo(m,i)
      end do

      do m=1,nO
        do e=1,nV
          r1(i,a) = r1(i,a) + t2(i,m,a,e)*cFov(m,e)
        end do
      end do

      do n=1,nO
        do f=1,nV
          r1(i,a) = r1(i,a) - t1(n,f)*OVOV(n,a,i,f)
        end do
      end do

      do m=1,nO
        do e=1,nV
          do f=1,nV
            r1(i,a) = r1(i,a) - 0.5d0*t2(i,m,e,f)*OVVV(m,a,e,f)
          end do
        end do
      end do

      do m=1,nO
        do n=1,nO
          do e=1,nV
            r1(i,a) = r1(i,a) - 0.5d0*t2(m,n,a,e)*OOVO(n,m,e,i)
          end do
        end do
      end do

    end do
  end do

! Final expression for t1 residue

  r1(:,:) = delta_ov(:,:)*t1(:,:) - r1(:,:)

end subroutine form_r1_nc
