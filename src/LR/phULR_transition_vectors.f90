subroutine phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

! Print transition vectors for linear response calculation

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
  double precision              :: dipole_int_aa(nBas,nBas,ncart)
  double precision              :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: c(nBas,nBas,nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Om(nSt)
  double precision,intent(in)   :: XpY(nSt,nSt)
  double precision,intent(in)   :: XmY(nSt,nSt)

! Local variables

  integer                       :: ia,jb,j,b
  integer                       :: maxS = 10
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: os(:)
  double precision,allocatable  :: S2(:)

! Memory allocation

  maxS = min(nSt,maxS)
  allocate(X(nSt),Y(nSt),os(maxS),S2(maxS))

! Compute oscillator strengths

  os(:) = 0d0
  if(ispin == 1) call phULR_oscillator_strength(nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,maxS, & 
                                                dipole_int_aa,dipole_int_bb,Om,XpY,XmY,os)

! Compute <S**2>

  call S2_expval(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,maxS,c,S,Om,XpY,XmY,S2)

! Print details about spin-conserved excitations

  if(ispin == 1) then

    do ia=1,maxS

      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

      print*,'-------------------------------------------------------------'
      write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A11,F6.4)') & 
              ' Excitation n. ',ia,': ',Om(ia)*HaToeV,' eV','  f = ',os(ia),'  <S**2> = ',S2(ia)
      print*,'-------------------------------------------------------------'

      ! Spin-up transitions

      jb = 0
      do j=nC(1)+1,nO(1)
        do b=nO(1)+1,nBas-nR(1)
          jb = jb + 1
          if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A -> ',b,'A = ',X(jb)
        end do
      end do
   
      jb = 0
      do j=nC(1)+1,nO(1)
        do b=nO(1)+1,nBas-nR(1)
          jb = jb + 1
          if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A <- ',b,'A = ',Y(jb)
        end do
      end do

      ! Spin-down transitions

      jb = 0
      do j=nC(2)+1,nO(2)
        do b=nO(2)+1,nBas-nR(2)
          jb = jb + 1
          if(abs(X(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'B -> ',b,'B = ',X(nSa+jb)
        end do
      end do
   
      jb = 0
      do j=nC(2)+1,nO(2)
        do b=nO(2)+1,nBas-nR(2)
          jb = jb + 1
          if(abs(Y(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'B <- ',b,'B = ',Y(nSa+jb)
        end do
      end do
     write(*,*)

    end do

  end if

! Print details about spin-flip excitations

  if(ispin == 2) then

    do ia=1,maxS

      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))


      print*,'-------------------------------------------------------------'
      write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A11,F6.4)') & 
              ' Excitation n. ',ia,': ',Om(ia)*HaToeV,' eV','  f = ',os(ia),'  <S**2> = ',S2(ia)
      print*,'-------------------------------------------------------------'

      ! Spin-up transitions

      jb = 0
      do j=nC(1)+1,nO(1)
        do b=nO(2)+1,nBas-nR(2)
          jb = jb + 1
          if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A -> ',b,'B = ',X(jb)
        end do
      end do
   
      jb = 0
      do j=nC(1)+1,nO(1)
        do b=nO(2)+1,nBas-nR(2)
          jb = jb + 1
          if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A <- ',b,'B = ',Y(jb)
        end do
      end do

      ! Spin-down transitions

      jb = 0
      do j=nC(2)+1,nO(2)
        do b=nO(1)+1,nBas-nR(1)
          jb = jb + 1
          if(abs(X(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A -> ',b,'B = ',X(nSa+jb)
        end do
      end do
   
      jb = 0
      do j=nC(2)+1,nO(2)
        do b=nO(1)+1,nBas-nR(1)
          jb = jb + 1
          if(abs(Y(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F10.6)') j,'A <- ',b,'B = ',Y(nSa+jb)
        end do
      end do
     write(*,*)

    end do

  end if

! Thomas-Reiche-Kuhn sum rule

  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
  write(*,*)

end subroutine 
