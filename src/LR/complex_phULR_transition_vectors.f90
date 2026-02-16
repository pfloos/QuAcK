subroutine complex_phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,&
                nCVS,nFC,occupations,virtuals,dipole_int_aa,dipole_int_bb,c,S,Om,XpY,XmY)

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
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas - minval(nO),nspin)
  complex*16                    :: dipole_int_aa(nBas,nBas,ncart)
  complex*16                    :: dipole_int_bb(nBas,nBas,ncart)
  complex*16,intent(in)         :: c(nBas,nBas,nspin)
  complex*16,intent(in)         :: S(nBas,nBas)
  complex*16,intent(in)         :: Om(nSt)
  complex*16,intent(in)         :: XpY(nSt,nSt)
  complex*16,intent(in)         :: XmY(nSt,nSt)

! Local variables

  integer                       :: ia,jb,j,b
  integer                       :: maxS = 30
  double precision,parameter    :: thres_vec = 0.1d0
  complex*16,allocatable        :: X(:)
  complex*16,allocatable        :: Y(:)
  complex*16,allocatable        :: os(:)
  complex*16,allocatable        :: S2(:)

! Memory allocation
  maxS = min(nSt,maxS)
  allocate(X(nSt),Y(nSt),os(maxS),S2(maxS))

! Not implemented for complex yet
!! Compute oscillator strengths
!
!  os(:) = 0d0
!  if(ispin == 1) call phULR_oscillator_strength(nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,maxS, & 
!                                                dipole_int_aa,dipole_int_bb,Om,XpY,XmY,os)
!
!! Compute <S**2>
!
!  call S2_expval(ispin,nBas,nC,nO,nV,nR,nS,nSa,nSb,nSt,maxS,c,S,Om,XpY,XmY,S2)

! Print details about spin-conserved excitations

  if(ispin == 1) then

    do ia=1,maxS

      X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
      Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

      print*,'-------------------------------------------------------------'
      write(*,'(A15,I3,A2,F15.6,A5,F15.6,A3)') & 
              ' Excitation n. ',ia,': ',real(Om(ia)*HaToeV),' + i ',aimag(Om(ia)*HaToeV),' eV'
      print*,'-------------------------------------------------------------'

      ! Spin-up transitions

      jb = 0
      do j=1,nO(1) - nFC(1)
        do b=nCVS(1)+1,nBas-nO(1)
          jb = jb + 1
          if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,1),'A -> ',&
            virtuals(b,1),'A = ',real(X(jb)),' + i ', aimag(X(jb))
        end do
      end do
   
      jb = 0
      do j=1,nO(1) - nFC(1)
        do b=nCVS(1)+1,nBas-nO(1)
          jb = jb + 1
          if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,1),'A <- ',&
               virtuals(b,1),'A = ',real(Y(jb)),' + i ', aimag(Y(jb))
        end do
      end do

      ! Spin-down transitions

      jb = 0
      do j=1,nO(2) - nFC(2)
        do b=nCVS(2)+1,nBas-nO(2)
          jb = jb + 1
          if(abs(X(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,2),'B -> ',&
              virtuals(b,2),'B = ',real(X(nSa+jb)),'+ i ',aimag(X(nSa + jb))
        end do
      end do
   
      jb = 0
      do j=1,nO(2)-nFC(2)
        do b=nCVS(2)+1,nBas-nO(2)
          jb = jb + 1
          if(abs(Y(nSa+jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,2),'B <- ',&
             virtuals(b,2),'B = ',real(Y(nSa+jb)),' + i ',aimag(Y(nSa + jb))
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


!      print*,'-------------------------------------------------------------'
!      write(*,'(A15,I3,A2,F10.6,A3,A6,F6.4,A11,F6.4)') & 
!              ' Excitation n. ',ia,': ',Om(ia)*HaToeV,' eV','  f = ',os(ia),'  <S**2> = ',S2(ia)
!      print*,'-------------------------------------------------------------'
 
      print*,'-------------------------------------------------------------'
      write(*,'(A15,I3,A2,F15.6,A5,F15.6,A3)') & 
              ' Excitation n. ',ia,': ',real(Om(ia)*HaToeV),' + i ',aimag(Om(ia)*HaToeV),' eV'
      print*,'-------------------------------------------------------------'

      ! Spin-up transitions

      jb = 0
      do j=1,nO(1) - nFC(1)
        do b=nCVS(2)+1,nBas-nO(2)
          jb = jb + 1
          if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,1),'A -> ',virtuals(b,2),&
          'B = ',real(X(jb)),' + i ', aimag(X(jb))
        end do
      end do
   
      jb = 0
      do j=1,nO(1) - nFC(1)
        do b=nCVS(2)+1,nBas-nO(2)
          jb = jb + 1
          if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,1),'A <- ',virtuals(b,2),&
          'B = ',real(Y(jb)),' + i ', aimag(Y(jb))
        end do
      end do

      ! Spin-down transitions

      jb = 0
      do j=1,nO(2) - nFC(2)
        do b=nCVS(1)+1,nBas-nO(1)
          jb = jb + 1
          if(abs(X(nSa + jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,2),'A -> ',virtuals(b,1),&
          'B = ',real(Y(jb)),' + i ', aimag(Y(jb))
        end do
      end do
   
      jb = 0
      do j=1,nO(2)-nFC(2)
        do b=nCVS(1)+1,nBas-nO(1)
          jb = jb + 1
          if(abs(Y(nSa + jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j,2),'A <- ',virtuals(b,1),&
          'B = ',real(Y(nSa+jb)),' + i ', aimag(Y(nSa+jb))
        end do
      end do
     write(*,*)

    end do

  end if

! Thomas-Reiche-Kuhn sum rule
!
!  write(*,'(A30,F10.6)') 'Thomas-Reiche-Kuhn sum rule = ',sum(os(:))
!  write(*,*)

end subroutine 
