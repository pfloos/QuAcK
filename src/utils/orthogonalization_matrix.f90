subroutine orthogonalization_matrix(ortho_type,nBas,S,X)

! Compute the orthogonalization matrix X

  implicit none

! Input variables

  integer,intent(in)            :: nBas,ortho_type
  double precision,intent(in)   :: S(nBas,nBas)

! Local variables

  logical                       :: debug
  double precision,allocatable  :: UVec(:,:),Uval(:)
  double precision,parameter    :: thresh = 1d-6

  integer                       :: i

! Output variables

  double precision,intent(out)  :: X(nBas,nBas)

  debug = .false.

! Type of orthogonalization ortho_type
!
!  1 = Lowdin
!  2 = Canonical
!  3 = SVD
!

  allocate(Uvec(nBas,nBas),Uval(nBas))

  if(ortho_type == 1) then

!   write(*,*)
!   write(*,*) ' Lowdin orthogonalization'
!   write(*,*)

    Uvec = S
    call diagonalize_matrix(nBas,Uvec,Uval)

    do i=1,nBas

      if(Uval(i) < thresh) then 

        write(*,*) 'Eigenvalue',i,' is very small in Lowdin orthogonalization = ',Uval(i)

      endif

      Uval(i) = 1d0/sqrt(Uval(i))

    enddo
    
    call ADAt(nBas,Uvec,Uval,X)

  elseif(ortho_type == 2) then

!   write(*,*)
!   write(*,*) 'Canonical orthogonalization'
!   write(*,*)

    Uvec = S
    call diagonalize_matrix(nBas,Uvec,Uval)

    do i=1,nBas

      if(Uval(i) > thresh) then 

        Uval(i) = 1d0/sqrt(Uval(i))

      else

        write(*,*) ' Eigenvalue',i,'too small for canonical orthogonalization'

      endif

    enddo
    
    call AD(nBas,Uvec,Uval)
    X = Uvec
 
  elseif(ortho_type == 3) then

!   write(*,*)
!   write(*,*) ' SVD-based orthogonalization NYI'
!   write(*,*)
   
!   Uvec = S
!   call diagonalize_matrix(nBas,Uvec,Uval)

!   do i=1,nBas
!     if(Uval(i) > thresh) then 
!       Uval(i) = 1d0/sqrt(Uval(i))
!     else
!       write(*,*) 'Eigenvalue',i,'too small for canonical orthogonalization'
!     endif
!   enddo
!   
!   call AD(nBas,Uvec,Uval)
!   X = Uvec
 
  endif

! Print results

  if(debug) then

    write(*,'(A28)') '----------------------'
    write(*,'(A28)') 'Orthogonalization matrix'
    write(*,'(A28)') '----------------------'
    call matout(nBas,nBas,X)
    write(*,*)

  endif

end subroutine orthogonalization_matrix
