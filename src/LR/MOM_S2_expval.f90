subroutine MOM_S2_expval(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,nCVS,nFC,occupations,virtuals,maxS,c,S,Omega,XpY,XmY,S2)

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
  integer,intent(in)            :: nCVS(nspin),nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
  integer,intent(in)            :: maxS
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Omega(nSt)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  integer                       :: m
  integer                       :: ia,i,a
  integer                       :: nOFC(nspin), nVCVS(nspin)
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
  
  nOFC(:)  = nO(:) - nFC(:)  ! occupied orbitals without counting frozen core
  nVCVS(:) = nV(:) - nCVS(:) ! Virtuals without counting CVS excluded ones

  allocate(OO(nOFC(1),nOFC(2)), OV(nOFC(1),nVCVS(2)),   &
           VO(nVCVS(1),nOFC(2)), VV(nVCVS(1),nVCVS(2)), &
           OOt(nOFC(2),nOFC(1)), OVt(nVCVS(2),nOFC(1)), &
           VOt(nOFC(2),nVCVS(1)), VVt(nVCVS(2),nVCVS(1)))

! Overlap matrix between spin-up and spin-down orbitals

  OO(:,:) = matmul(transpose(c(:,occupations(1: nOFC(1),1),1)),  &
            matmul(S,c(:,occupations(1: nOFC(2),2),2)))
  OV(:,:) = matmul(transpose(c(:,occupations(1: nOFC(1),1),1)),  &
            matmul(S,c(:,virtuals(1:nVCVS(2),2),2)))
  VO(:,:) = matmul(transpose(c(:,virtuals(1:nVCVS(1),1),1)),   &
            matmul(S,c(:,occupations(1: nOFC(2),2),2)))
  VV(:,:) = matmul(transpose(c(:,virtuals(1:nVCVS(1),1),1)),  &
            matmul(S,c(:,virtuals(1:nVCVS(2),2),2)))

  OOt(:,:) = transpose(OO(:,:))
  OVt(:,:) = transpose(OV(:,:))
  VOt(:,:) = transpose(VO(:,:))
  VVt(:,:) = transpose(VV(:,:))

!-------------------------!
! <S**2> for ground state !
!-------------------------!

  S2_exact = dble(nO(1) - nO(2))/2d0*(dble(nO(1) - nO(2))/2d0 + 1d0)
  S2_gs    = S2_exact + dble(nO(2)) - sum(OO(:,:)**2)

!------------------------------------------!
! <S**2> for spin-conserved-excited states !
!------------------------------------------!

  if(ispin == 1) then

    allocate(Xa (nOFC(1) ,nVCVS(1)), Ya (nOFC(1) ,nVCVS(1)),  &
             Xb (nOFC(2) ,nVCVS(2)), Yb (nOFC(2) ,nVCVS(2)),  &
             Xat(nVCVS(1),nOFC(1)) , Yat(nVCVS(1),nOFC(1)),   &
             Xbt(nVCVS(2),nOFC(2)) , Ybt(nVCVS(2),nOFC(2)))
    
    do m=1,maxS

      ia = 0
      do i=1,nOFC(1)
        do a=1,nVCVS(1)
          ia = ia + 1
          Xa(i,a) = 0.5d0*(XpY(m,ia) + XmY(m,ia))
          Ya(i,a) = 0.5d0*(XpY(m,ia) - XmY(m,ia))
        end do
      end do

      ia = 0
      do i=1,nOFC(2)
        do a=1,nVCVS(2)
          ia = ia + 1
          Xb(i,a) = 0.5d0*(XpY(m,nSa+ia) + XmY(m,nSa+ia))
          Yb(i,a) = 0.5d0*(XpY(m,nSa+ia) - XmY(m,nSa+ia))
        end do
      end do

      Xat(:,:) = transpose(Xa(:,:))
      Xbt(:,:) = transpose(Xb(:,:))
      Yat(:,:) = transpose(Ya(:,:))
      Ybt(:,:) = transpose(Yb(:,:))

      S2(m) = S2_gs                                                             &                                                         
            + trace_matrix(nVCVS(1),matmul(Xat,matmul(OO, matmul(OOt,Xa))))     &   
            + trace_matrix(nVCVS(2),matmul(Xbt,matmul(OOt,matmul(OO,Xb))))      &
            - trace_matrix(nOFC(1), matmul(Xa, matmul(VO, matmul(VOt,Xat))))    &
            - trace_matrix(nOFC(2), matmul(Xb, matmul(OVt,matmul(OV,Xbt))))     &
      - 2d0 * trace_matrix(nOFC(1), matmul(OO, matmul(Xb, matmul(VVt,Xat))))    &
            
      - 2d0 * trace_matrix(nVCVS(2),matmul(OVt,matmul(Xa, matmul(VO,Yb))))      &
      - 2d0 * trace_matrix(nVCVS(1),matmul(VO,matmul(Xb,  matmul(OVt,Ya))))     &
            
            - trace_matrix(nVCVS(1),matmul(Yat,matmul(OO, matmul(OOt,Ya))))     &
            - trace_matrix(nVCVS(2),matmul(Ybt,matmul(OOt,matmul(OO,Yb))))      &
            + trace_matrix(nOFC(1), matmul(Ya,matmul(VO,  matmul(VOt,Yat))))    &
            + trace_matrix(nOFC(2), matmul(Yb,matmul(OVt, matmul(OV,Ybt))))     &
      + 2d0 * trace_matrix(nOFC(1), matmul(Ya,matmul(VV,  matmul(Ybt,OOt))))

    end do

  end if

!------------------------------------------!
! <S**2> for spin-conserved-excited states !
!------------------------------------------!

  if(ispin == 2) then

    allocate(Xa(nOFC(1),nVCVS(2)), Ya(nOFC(1),nVCVS(2)),    &
             Xb(nOFC(2),nVCVS(1)), Yb(nOFC(2),nVCVS(1)),    &
             Xat(nVCVS(2),nOFC(1)),Yat(nVCVS(2),nOFC(1)),   &
             Xbt(nVCVS(1),nOFC(2)),Ybt(nVCVS(1),nOFC(2)))

    do m=1,maxS

      ia = 0
      do i=1,nOFC(1)
        do a=1,nVCVS(2)
          ia = ia + 1
          Xa(i,a) = 0.5d0*(XpY(m,ia) + XmY(m,ia))
          Ya(i,a) = 0.5d0*(XpY(m,ia) - XmY(m,ia))
        end do
      end do

      ia = 0
      do i=1,nOFC(2)
        do a=1,nVCVS(1)
          ia = ia + 1
          Xb(i,a) = 0.5d0*(XpY(m,nSa+ia) + XmY(m,nSa+ia))
          Yb(i,a) = 0.5d0*(XpY(m,nSa+ia) - XmY(m,nSa+ia))
        end do
      end do

      Xat(:,:) = transpose(Xa(:,:))
      Xbt(:,:) = transpose(Xb(:,:))
      Yat(:,:) = transpose(Ya(:,:))
      Ybt(:,:) = transpose(Yb(:,:))

      S2(m) = S2_gs + dble(nO(2) - nO(1)) + 1d0  
        
      S2(m) = S2(m)                                                              &
            + trace_matrix(nVCVS(1),matmul(Xbt,matmul(OOt,matmul(OO,Xb))))       &     
            - trace_matrix(nOFC(2), matmul(Xb, matmul(VO, matmul(VOt,Xbt))))     &     
            + trace_matrix(nOFC(2), matmul(Xb,VO))**2                            &     
            + trace_matrix(nVCVS(2),matmul(Yat,matmul(OO, matmul(OOt,Ya))))      &     
            + trace_matrix(nOFC(1), matmul(Ya, matmul(OVt, matmul(OV,Yat))))     &     
            + trace_matrix(nOFC(1), matmul(Ya,OVt))**2                           &      
      - 2d0 * trace_matrix(nOFC(2), matmul(Xb,VO))                               &                              
            * trace_matrix(nOFC(1), matmul(Ya,OVt))                              &

            + trace_matrix(nVCVS(2),matmul(Xat,matmul(OO,matmul(OOt,Xa))))       &      
            - trace_matrix(nOFC(1),matmul(Xa,matmul(OVt,matmul(OV,Xat))))        &      
            + trace_matrix(nOFC(1),matmul(Xa,OVt))**2                            &      
            + trace_matrix(nVCVS(1),matmul(Ybt,matmul(OOt,matmul(OO,Yb))))       &      
            - trace_matrix(nOFC(2),matmul(Yb,matmul(VO,matmul(VOt,Ybt))))        &      
            + trace_matrix(nVCVS(1),matmul(Ybt,VOt))**2                          &      
      - 2d0 * trace_matrix(nOFC(1),matmul(Xa,OVt))                               &      
            * trace_matrix(nOFC(2),matmul(Yb,VO))

    end do

  end if

end subroutine 
