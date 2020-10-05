subroutine unrestricted_S2_expval(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,maxS,c,S,XpY,XmY,S2)

! Compute <S**2> for linear response excited states

  implicit none
  include 'parameters.h'

! Input variables


  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  integer,intent(in)            :: maxS
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  integer                       :: m
  integer                       :: ia,i,a
  double precision              :: S2_exact
  double precision              :: S2_gs
  double precision,allocatable  :: Xa(:,:), Xb(:,:), Ya(:,:), Yb(:,:)
  double precision,allocatable  :: Xat(:,:),Xbt(:,:),Yat(:,:),Ybt(:,:)
  double precision,allocatable  :: OO(:,:), OV(:,:), VO(:,:), VV(:,:)
  double precision,allocatable  :: OOt(:,:),OVt(:,:),VOt(:,:),VVt(:,:)
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: S2(maxS)

! Memory allocation

  allocate(OO(nO(1)-nC(1),nO(2)-nC(2)), OV(nO(1)-nC(1),nV(2)-nR(2)), VO(nV(1)-nR(1),nO(2)-nC(2)), VV(nV(1)-nR(1),nV(2)-nR(2)), &
           OOt(nO(2)-nC(2),nO(1)-nC(1)),OVt(nV(2)-nR(2),nO(1)-nC(1)),VOt(nO(2)-nC(2),nV(1)-nR(1)),VVt(nV(2)-nR(2),nV(1)-nR(1)))

! Overlap matrix between spin-up and spin-down orbitals

  OO(:,:) = matmul(transpose(c(:,nC(1)+1:nO(1)     ,1)),matmul(S,c(:,nC(2)+1:nO(2)     ,2)))
  OV(:,:) = matmul(transpose(c(:,nC(1)+1:nO(1)     ,1)),matmul(S,c(:,nO(2)+1:nBas-nR(2),2)))
  VO(:,:) = matmul(transpose(c(:,nO(1)+1:nBas-nR(1),1)),matmul(S,c(:,nC(2)+1:nO(2)     ,2)))
  VV(:,:) = matmul(transpose(c(:,nO(1)+1:nBas-nR(1),1)),matmul(S,c(:,nO(2)+1:nBas-nR(2),2)))

  OOt(:,:) = transpose(OO(:,:))
  OVt(:,:) = transpose(OV(:,:))
  VOt(:,:) = transpose(VO(:,:))
  VVt(:,:) = transpose(VV(:,:))

!-------------------------!
! <S**2> for ground state !
!-------------------------!

  S2_exact = dble(nO(1) - nO(2))/2d0*(dble(nO(1) - nO(2))/2d0 + 1d0)
  S2_gs    = S2_exact + nO(2) - sum(OO(:,:)**2)

!------------------------------------------!
! <S**2> for spin-conserved-excited states !
!------------------------------------------!

  if(ispin == 1) then

    allocate(Xa(nO(1)-nC(1),nV(1)-nR(1)), Ya(nO(1)-nC(1),nV(1)-nR(1)), Xb(nO(2)-nC(2),nV(2)-nR(2)), Yb(nO(2)-nC(2),nV(2)-nR(2)), &
             Xat(nV(1)-nR(1),nO(1)-nC(1)),Yat(nV(1)-nR(1),nO(1)-nC(1)),Xbt(nV(2)-nR(2),nO(2)-nC(2)),Ybt(nV(2)-nR(2),nO(2)-nC(2)))

    do m=1,maxS

      ia = 0
      do i=nC(1)+1,nO(1)
        do a=1,nV(1)-nR(1)
          ia = ia + 1
          Xa(i,a) = 0.5d0*(XpY(m,ia) + XmY(m,ia))
          Ya(i,a) = 0.5d0*(XpY(m,ia) - XmY(m,ia))
        end do
      end do

      ia = 0
      do i=nC(2)+1,nO(2)
        do a=1,nV(2)-nR(2)
          ia = ia + 1
          Xb(i,a) = 0.5d0*(XpY(m,nSa+ia) + XmY(m,nSa+ia))
          Yb(i,a) = 0.5d0*(XpY(m,nSa+ia) - XmY(m,nSa+ia))
        end do
      end do

      Xat(:,:) = transpose(Xa(:,:))
      Xbt(:,:) = transpose(Xb(:,:))
      Yat(:,:) = transpose(Ya(:,:))
      Ybt(:,:) = transpose(Yb(:,:))

      S2(m) = S2_gs &
            + trace_matrix(nV(1),matmul(Xat,matmul(OO,matmul(OOt,Xa)))) &
            + trace_matrix(nV(2),matmul(Xbt,matmul(OOt,matmul(OO,Xb)))) &
            - trace_matrix(nO(1),matmul(Xa,matmul(VO,matmul(VOt,Xat)))) &
            - trace_matrix(nO(2),matmul(Xb,matmul(OVt,matmul(OV,Xbt)))) &
            - 2d0*trace_matrix(nO(1),matmul(OO,matmul(Xb,matmul(VVt,Xat)))) &
            
            - 2d0*trace_matrix(nV(2),matmul(OVt,matmul(Xa,matmul(VO,Yb)))) &
            - 2d0*trace_matrix(nV(1),matmul(VO,matmul(Xb,matmul(OVt,Ya)))) &
            
            - trace_matrix(nV(1),matmul(Yat,matmul(OO,matmul(OOt,Ya)))) &
            - trace_matrix(nV(2),matmul(Ybt,matmul(OOt,matmul(OO,Yb)))) &
            + trace_matrix(nO(1),matmul(Ya,matmul(VO,matmul(VOt,Yat)))) &
            + trace_matrix(nO(2),matmul(Yb,matmul(OVt,matmul(OV,Ybt)))) &
            + 2d0*trace_matrix(nO(1),matmul(Ya,matmul(VV,matmul(Ybt,OOt))))

    end do

  end if

!------------------------------------------!
! <S**2> for spin-conserved-excited states !
!------------------------------------------!

  if(ispin == 2) then

    allocate(Xa(nO(1)-nC(1),nV(2)-nR(2)), Ya(nO(1)-nC(1),nV(2)-nR(2)), Xb(nO(2)-nC(2),nV(1)-nR(1)), Yb(nO(2)-nC(2),nV(1)-nR(1)), &
             Xat(nV(2)-nR(2),nO(1)-nC(1)),Yat(nV(2)-nR(2),nO(1)-nC(1)),Xbt(nV(1)-nR(1),nO(2)-nC(2)),Ybt(nV(1)-nR(1),nO(2)-nC(2)))

    do m=1,maxS

      ia = 0
      do i=nC(1)+1,nO(1)
        do a=1,nV(2)-nR(2)
          ia = ia + 1
          Xa(i,a) = 0.5d0*(XpY(m,ia) + XmY(m,ia))
          Ya(i,a) = 0.5d0*(XpY(m,ia) - XmY(m,ia))
        end do
      end do

      ia = 0
      do i=nC(2)+1,nO(2)
        do a=1,nV(1)-nR(1)
          ia = ia + 1
          Xb(i,a) = 0.5d0*(XpY(m,nSa+ia) + XmY(m,nSa+ia))
          Yb(i,a) = 0.5d0*(XpY(m,nSa+ia) - XmY(m,nSa+ia))
        end do
      end do

      Xat(:,:) = transpose(Xa(:,:))
      Xbt(:,:) = transpose(Xb(:,:))
      Yat(:,:) = transpose(Ya(:,:))
      Ybt(:,:) = transpose(Yb(:,:))

      S2(m) = S2_gs & 

            + trace_matrix(nV(1),matmul(Xbt,matmul(OOt,matmul(OO,Xb)))) &
            - trace_matrix(nO(2),matmul(Xb,matmul(VO,matmul(VOt,Xbt)))) &
            + trace_matrix(nO(2),matmul(Xb,VO))**2 &
            + trace_matrix(nV(2),matmul(Yat,matmul(OO,matmul(OOt,Ya)))) &
            + trace_matrix(nO(1),matmul(Ya,matmul(OVt,matmul(OV,Yat)))) &
            + trace_matrix(nO(1),matmul(Ya,OVt))**2 &
            - 2d0*trace_matrix(nO(2),matmul(Xb,VO))*trace_matrix(nO(1),matmul(Ya,OVt)) &

            + trace_matrix(nV(2),matmul(Xat,matmul(OO,matmul(OOt,Xa)))) &
            - trace_matrix(nO(1),matmul(Xa,matmul(OVt,matmul(OV,Xat)))) &
            + trace_matrix(nO(1),matmul(Xa,OVt))**2 &
            + trace_matrix(nV(1),matmul(Ybt,matmul(OOt,matmul(OO,Yb)))) &
            - trace_matrix(nO(2),matmul(Yb,matmul(VO,matmul(VOt,Ybt)))) &
            + trace_matrix(nV(1),matmul(Ybt,VOt))**2 &
            - 2d0*trace_matrix(nO(1),matmul(Xa,OVt))*trace_matrix(nO(2),matmul(Yb,VO))

    end do

  end if

end subroutine unrestricted_S2_expval
