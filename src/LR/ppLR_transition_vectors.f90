subroutine ppLR_transition_vectors(spin_allowed,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

! Print transition vectors for p-p linear response calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: spin_allowed
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: X1(nVV,nVV)
  double precision,intent(in)   :: Y1(nOO,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: X2(nVV,nOO)
  double precision,intent(in)   :: Y2(nOO,nOO)

! Local variables

  integer                       :: a,b,c,d,ab,cd
  integer                       :: i,j,k,l,ij,kl
  integer                       :: maxOO = 25
  integer                       :: maxVV = 25
  double precision              :: S2
  double precision,parameter    :: thres_vec = 0.1d0
  double precision,allocatable  :: os1(:)
  double precision,allocatable  :: os2(:)

! Memory allocation

  maxOO = min(nOO,maxOO)
  maxVV = min(nVV,maxVV)

  allocate(os1(nVV),os2(nOO))

! Compute oscillator strengths

  os1(:) = 0d0
  os2(:) = 0d0

  if(spin_allowed) & 
    call ppLR_oscillator_strength(nBas,nC,nO,nV,nR,nOO,nVV,maxOO,maxVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2,os1,os2)

!-----------------------------------------------!
! Print details about excitations for pp sector !
!-----------------------------------------------!

  write(*,*) '*****************************'
  write(*,*) '*** (N+2)-electron states ***'
  write(*,*) '*****************************'
  write(*,*)

  do ab=1,maxVV

    ! <S**2> values

    if(spin_allowed) then 
      S2 = 0d0
    else
      S2 = 2d0
    end if

    print*,'-------------------------------------------------------------'
    write(*,'(A20,I3,A2,F12.6,A3,A6,F6.4,A11,F6.4)') &
            ' p-p excitation n. ',ab,': ',Om1(ab)*HaToeV,' eV','  f = ',os1(ab),'  <S**2> = ',S2
    print*,'-------------------------------------------------------------'

   if(spin_allowed) then

      cd = 0
      do c=nO+1,nBas-nR
        do d=c,nBas-nR
          cd = cd + 1 
          if(abs(X1(cd,ab)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') c,' -- ',d,' = ',X1(cd,ab)/sqrt(2d0)
        end do
      end do
     
      kl = 0
      do k=nC+1,nO
        do l=k,nO
          kl = kl + 1
          if(abs(Y1(kl,ab)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') k,' -- ',l,' = ',Y1(kl,ab)/sqrt(2d0)
        end do
      end do

    else

      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1 
          if(abs(X1(cd,ab)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') c,' -- ',d,' = ',X1(cd,ab)/sqrt(2d0)
        end do
      end do
     
      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1
          if(abs(Y1(kl,ab)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') k,' -- ',l,' = ',Y1(kl,ab)/sqrt(2d0)
        end do
      end do

    end if

   write(*,*)

  end do

! Thomas-Reiche-Kuhn sum rule

!  if(nVV > 0) write(*,'(A50,F10.6)') 'Thomas-Reiche-Kuhn sum rule for p-p sector = ',sum(os1(:))
!  write(*,*)

!-----------------------------------------------!
! Print details about excitations for hh sector !
!-----------------------------------------------!

  write(*,*) '*****************************'
  write(*,*) '*** (N-2)-electron states ***'
  write(*,*) '*****************************'
  write(*,*)

  do ij=nOO,nOO-maxOO+1,-1

    ! <S**2> values

    if(spin_allowed) then 
      S2 = 0d0
    else
      S2 = 2d0
    end if

    print*,'-------------------------------------------------------------'
    write(*,'(A20,I3,A2,F12.6,A3,A6,F6.4,A11,F6.4)') &
            ' h-h excitation n. ',ij,': ',Om2(ij)*HaToeV,' eV','  f = ',os2(ij),'  <S**2> = ',S2
    print*,'-------------------------------------------------------------'

   if(spin_allowed) then

      cd = 0
      do c=nO+1,nBas-nR
        do d=c,nBas-nR
          cd = cd + 1
          if(abs(X2(cd,ij)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') c,' -- ',d,' = ',X2(cd,ij)/sqrt(2d0)
        end do
      end do
     
      kl = 0
      do k=nC+1,nO
        do l=k,nO
          kl = kl + 1
          if(abs(Y2(kl,ij)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') k,' -- ',l,' = ',Y2(kl,ij)/sqrt(2d0)
        end do
      end do

    else

      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1
          if(abs(X2(cd,ij)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') c,' -- ',d,' = ',X2(cd,ij)/sqrt(2d0)
        end do
      end do
     
      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1
          if(abs(Y2(kl,ij)) > thres_vec) write(*,'(I3,A4,I3,A3,F10.6)') k,' -- ',l,' = ',Y2(kl,ij)/sqrt(2d0)
        end do
      end do

    end if

   write(*,*)

  end do

! Thomas-Reiche-Kuhn sum rule

!  if(nOO > 0) write(*,'(A50,F10.6)') 'Thomas-Reiche-Kuhn sum rule for h-h sector = ',sum(os2(:))
!  write(*,*)

end subroutine 
