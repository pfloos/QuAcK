subroutine form_FB(nC,nO,nV,nR,foo,fvv,fov,OOOV,OOVV,OVVV,t,FooB,FvvB,FovB)

! Form the effective Fock operator for BCCD

  implicit none

! Input variables

  integer,intent(in)            :: nC,nO,nV,nR
  double precision,intent(in)   :: foo(nO-nC,nO-nC)
  double precision,intent(in)   :: fvv(nV-nR,nV-nR)
  double precision,intent(in)   :: fov(nO-nC,nV-nR)
  double precision,intent(in)   :: OOOV(nO-nC,nO-nC,nO-nC,nV-nR)
  double precision,intent(in)   :: OOVV(nO-nC,nO-nC,nV-nR,nV-nR)
  double precision,intent(in)   :: OVVV(nO-nC,nV-nR,nV-nR,nV-nR)
  double precision,intent(in)   :: t(nO-nC,nO-nC,nV-nR,nV-nR)

! Local variables

  integer                       :: i,j,k,a,b,c

! Output variables

  double precision,intent(out)  :: FooB(nO-nC,nO-nC)
  double precision,intent(out)  :: FvvB(nV-nR,nV-nR)
  double precision,intent(out)  :: FovB(nO-nC,nV-nR)

! Occupied-occupied block

  FooB(:,:) = foo(:,:)
  do i=1,nO-nC
    do j=1,nO-nC
      do k=1,nO-nC
        do a=1,nV-nR
          do b=1,nV-nR
            FooB(i,k) = FooB(i,k) + 0.5d0*OOVV(k,j,a,b)*t(i,j,a,b)
          enddo
        enddo
      enddo
    enddo
  enddo

! Virtual-virtual block  

  FvvB(:,:) = fvv(:,:)
  do a=1,nV-nR
    do b=1,nV-nR
      do i=1,nO-nC
        do j=1,nO-nC
          do c=1,nV-nR
            FvvB(a,c) = FvvB(a,c) - 0.5d0*OOVV(i,j,c,b)*t(i,j,a,b)
          enddo
        enddo
      enddo
    enddo
  enddo

! Occupied-virtual block  

  FovB(:,:) = fov(:,:)
  do i=1,nO-nC
    do a=1,nV-nR

      do j=1,nO-nC
        do b=1,nV-nR
          FovB(i,a) = FovB(i,a) - fov(j,b)*t(i,j,a,b)
        enddo
      enddo

      do j=1,nO-nC
        do b=1,nV-nR
          do c=1,nV-nR
            FovB(i,a) = FovB(i,a) + 0.5d0*OVVV(j,a,b,c)*t(i,j,b,c)
          enddo
        enddo
      enddo

      do j=1,nO-nC
        do k=1,nO-nC
          do b=1,nV-nR
            FovB(i,a) = FovB(i,a) - 0.5d0*OOOV(j,k,i,b)*t(j,k,a,b)
          enddo
        enddo
      enddo

    enddo
  enddo

end subroutine form_FB
