subroutine ppRLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

  ! 
  ! Solve the pp-RPA linear eigenvalue problem
  !
  ! right eigen-problem: H   R = R w
  ! left  eigen-problem: H.T L = L w
  !
  !     where L.T R = 1
  !
  !
  !        (+C    +B) 
  !   H =  (        )  where C = C.T and D = D.T
  !        (-B.T  -D)
  !
  !      (w1   0)        (X1   X2)          (+X1  +X2)
  !  w = (      ),   R = (       )  and L = (        )
  !      (0   w2)        (Y1   Y2)          (-Y1  -Y2)
  ! 
  !
  !  the normalisation condition reduces to
  !
  !    X1.T X2 - Y1.T Y2 = 0
  !    X1.T X1 - Y1.T Y1 = 1
  !    X2.T X2 - Y2.T Y2 = 1
  !

  implicit none
  include 'parameters.h'

  logical,intent(in)           :: TDA
  integer,intent(in)           :: nOO,nVV
  double precision,intent(in)  :: Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO)
  double precision,intent(out) :: Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV)
  double precision,intent(out) :: Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO)
  double precision,intent(out) :: EcRPA
  
  logical                       :: imp_bio,verbose
  integer                       :: i,j
  integer                       :: Npp
  double precision              :: EcRPA1,EcRPA2
  double precision              :: thr_d,thr_nd,thr_deg
  double precision,allocatable  :: M(:,:),Z(:,:),Om(:)
  
  double precision,allocatable  :: U(:,:),VT(:,:),S(:)
  double precision,allocatable  :: work(:)
  integer                       :: lwork, info

  double precision,external     :: trace_matrix

  nPP = nOO + nVV

  allocate(M(nPP,nPP),Z(nPP,nPP),Om(nPP))
  allocate(U(nPP,nPP),VT(nPP,nPP),S(nPP))

  ! Hermitian case for TDA

  if(TDA) then

    X1(:,:) = +Cpp(:,:)
    Y1(:,:) = 0d0
    if(nVV > 0) call diagonalize_matrix(nVV,X1,Om1)

    X2(:,:) = 0d0
    Y2(:,:) = -Dpp(:,:)
    if(nOO > 0) call diagonalize_matrix(nOO,Y2,Om2)

  else

    ! Diagonal blocks 

    M(    1:nVV,    1:nVV) = + Cpp(1:nVV,1:nVV)
    M(nVV+1:nPP,nVV+1:nPP) = - Dpp(1:nOO,1:nOO)

    ! Off-diagonal blocks

    M(    1:nVV,nVV+1:nPP) = +           Bpp(1:nVV,1:nOO)
    M(nVV+1:nPP,    1:nVV) = - transpose(Bpp(1:nVV,1:nOO))

    allocate(work(1))
    lwork = -1
    call dgesvd('A','A',nPP,nPP,M(1,1),nPP,S(1),U(1,1),nPP,VT(1,1),nPP,work,lwork,info)
    lwork = int(work(1))
    deallocate(work)
    allocate(work(lwork))
    
    call dgesvd('A','A',nPP,nPP,M(1,1),nPP,S(1),U(1,1),nPP,VT(1,1),nPP,work,lwork,info)

    write (*,*) info
    call vecout(nPP,S)

    ! Diagonal blocks 

    M(    1:nVV,    1:nVV) = + Cpp(1:nVV,1:nVV)
    M(nVV+1:nPP,nVV+1:nPP) = - Dpp(1:nOO,1:nOO)

    ! Off-diagonal blocks

    M(    1:nVV,nVV+1:nPP) = +           Bpp(1:nVV,1:nOO)
    M(nVV+1:nPP,    1:nVV) = - transpose(Bpp(1:nVV,1:nOO))

!   if((nOO .eq. 0) .or. (nVV .eq. 0)) then

      ! Diagonalize the pp matrix

      if(nPP > 0) call diagonalize_general_matrix(nPP,M,Om,Z)

      ! Split the various quantities in ee and hh parts

      call sort_ppRPA(nOO,nVV,nPP,Om,Z,Om1,X1,Y1,Om2,X2,Y2)

   ! else

   !     thr_d   = 1d-6   ! to determine if     diagonal elements of L.T x R are close enouph to 1
   !     thr_nd  = 1d-6   ! to determine if non-diagonal elements of L.T x R are close enouph to 1
   !     thr_deg = 1d-8   ! to determine if two eigenvectors are degenerate or not
   !     imp_bio = .True. ! impose bi-orthogonality
   !     verbose = .False.
   !     call diagonalize_nonsym_matrix(Npp, M, Z, Om, thr_d, thr_nd, thr_deg, imp_bio, verbose)
   
   !     do i = 1, nOO
   !       Om2(i) = Om(i)
   !       do j = 1, nVV
   !         X2(j,i) = Z(j,i)
   !       enddo
   !       do j = 1, nOO
   !         Y2(j,i) = Z(nVV+j,i)
   !       enddo
   !     enddo
   
   !     do i = 1, nVV
   !       Om1(i) = Om(nOO+i)
   !       do j = 1, nVV
   !         X1(j,i) = M(j,nOO+i)
   !       enddo
   !       do j = 1, nOO
   !         Y1(j,i) = M(nVV+j,nOO+i)
   !       enddo
   !     enddo

   ! endif

  end if

  ! Compute the RPA correlation energy

  EcRPA = 0.5d0 * (sum(Om1) - sum(Om2) - trace_matrix(nVV,Cpp) - trace_matrix(nOO,Dpp))
  EcRPA1 = +sum(Om1) - trace_matrix(nVV,Cpp)
  EcRPA2 = -sum(Om2) - trace_matrix(nOO,Dpp)

  ! if(abs(EcRPA - EcRPA1) > 1d-6 .or. abs(EcRPA - EcRPA2) > 1d-6) then
  !   print*,'!!! Issue in pp-RPA linear reponse calculation RPA1 != RPA2 !!!'
  ! endif

  deallocate(M,Z,Om)

end subroutine 
