subroutine mo_guess(nBas,guess_type,S,Hc,X,c)

!  Guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)

! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

  if(guess_type == 1) then

    call core_guess(nBas,Hc,X,c)

  elseif(guess_type == 2) then

    call huckel_guess(nBas,S,Hc,X,c)

  elseif(guess_type == 3) then

    call random_number(c)

    c(:,:) = 2d0*c(:,:) - 1d0

  else

    print*,'Wrong guess option'
    stop

  endif

end subroutine
