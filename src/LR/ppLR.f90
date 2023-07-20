subroutine ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

! Compute the pp channel of the linear response: see Scuseria et al. JCP 139, 104113 (2013)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: Bpp(nVV,nOO)
  double precision,intent(in)   :: Cpp(nVV,nVV)
  double precision,intent(in)   :: Dpp(nOO,nOO)
  
! Local variables

  double precision              :: trace_matrix
  double precision              :: EcRPA1
  double precision              :: EcRPA2
  double precision,allocatable  :: M(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: Om(:)

! Output variables

  double precision,intent(out)  :: Om1(nVV)
  double precision,intent(out)  :: X1(nVV,nVV)
  double precision,intent(out)  :: Y1(nOO,nVV)
  double precision,intent(out)  :: Om2(nOO)
  double precision,intent(out)  :: X2(nVV,nOO)
  double precision,intent(out)  :: Y2(nOO,nOO)
  double precision,intent(out)  :: EcRPA

! Memory allocation

  allocate(M(nOO+nVV,nOO+nVV),Z(nOO+nVV,nOO+nVV),Om(nOO+nVV))

!-------------------------------------------------!
! Solve the p-p eigenproblem                      !
!-------------------------------------------------!
!                                                 !
!  |  C   B | | X1  X2 |   | w1  0  | | X1  X2 |  !
!  |        | |        | = |        | |        |  !
!  | -Bt -D | | Y1  Y2 |   | 0   w2 | | Y1  Y2 |  !
!                                                 !
!-------------------------------------------------!

  if(TDA) then

    X1(:,:) = +Cpp(:,:)
    Y1(:,:) = 0d0
    if(nVV > 0) call diagonalize_matrix(nVV,X1,Om1)

    X2(:,:) = 0d0
    Y2(:,:) = -Dpp(:,:)
    if(nOO > 0) call diagonalize_matrix(nOO,Y2,Om2)

  else

  ! Diagonal blocks 

    M(    1:nVV    ,    1:nVV)     = + Cpp(1:nVV,1:nVV)
    M(nVV+1:nVV+nOO,nVV+1:nVV+nOO) = - Dpp(1:nOO,1:nOO)

  ! Off-diagonal blocks

    M(    1:nVV    ,nVV+1:nOO+nVV) = -           Bpp(1:nVV,1:nOO)
    M(nVV+1:nOO+nVV,    1:nVV)     = + transpose(Bpp(1:nVV,1:nOO))

  ! Diagonalize the p-p matrix

    if(nOO+nVV > 0) call diagonalize_general_matrix(nOO+nVV,M,Om,Z)

  ! Split the various quantities in p-p and h-h parts

    call sort_ppRPA(nOO,nVV,Om,Z,Om1,X1,Y1,Om2,X2,Y2)

  end if

! Compute the RPA correlation energy

  EcRPA = 0.5d0*( sum(Om1) - sum(Om2) - trace_matrix(nVV,Cpp) - trace_matrix(nOO,Dpp) )
  EcRPA1 = +sum(Om1) - trace_matrix(nVV,Cpp)
  EcRPA2 = -sum(Om2) - trace_matrix(nOO,Dpp)

  if(abs(EcRPA - EcRPA1) > 1d-6 .or. abs(EcRPA - EcRPA2) > 1d-6) & 
    print*,'!!! Issue in pp-RPA linear reponse calculation RPA1 != RPA2 !!!'

end subroutine 
