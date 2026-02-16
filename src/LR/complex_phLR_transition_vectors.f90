subroutine complex_phLR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,&
                nCVS,nFC,occupations,virtuals,Om,XpY,XmY)

! Print transition vectors for linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nCVS
  integer,intent(in)            :: nFC
  integer,intent(in)            :: occupations(nO-nFC)
  integer,intent(in)            :: virtuals(nBas - nO)
  complex*16,intent(in)         :: Om(nS)
  complex*16,intent(in)         :: XpY(nS,nS)
  complex*16,intent(in)         :: XmY(nS,nS)

! Local variables

  integer                       :: ia,jb,j,b
  integer                       :: maxS = 30
  double precision,parameter    :: thres_vec = 0.1d0
  complex*16,allocatable        :: X(:)
  complex*16,allocatable        :: Y(:)

! Memory allocation
  maxS = min(nS,maxS)
  allocate(X(nS),Y(nS))

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


  do ia=1,maxS

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    print*,'-------------------------------------------------------------'
    write(*,'(A15,I3,A2,F15.6,A5,F15.6,A3)') & 
            ' Excitation n. ',ia,': ',real(Om(ia)*HaToeV),' + i ',aimag(Om(ia)*HaToeV),' eV'
    print*,'-------------------------------------------------------------'

    jb = 0
    do j=1,nO - nFC
      do b=nCVS+1,nBas-nO
        jb = jb + 1
        if(abs(X(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j),' -> ',&
          virtuals(b),' = ',real(X(jb)),' + i ', aimag(X(jb))
      end do
    end do
 
    jb = 0
    do j=1,nO - nFC
      do b=nCVS+1,nBas-nO
        jb = jb + 1
        if(abs(Y(jb)) > thres_vec) write(*,'(I3,A5,I3,A4,F15.6,A5,F15.6)') occupations(j),' <- ',&
             virtuals(b),' = ',real(Y(jb)),' + i ', aimag(Y(jb))
      end do
    end do

  end do
end subroutine 
