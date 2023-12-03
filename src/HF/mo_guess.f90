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

  double precision,intent(inout)  :: c(nBas,nBas)

  if(guess_type == 0) then

    write(*,*) 'Reading HF coefficients...'
    c(:,:) = c(:,:)

  elseif(guess_type == 1) then

    write(*,*) 'Core guess...'
    call core_guess(nBas,Hc,X,c)

  elseif(guess_type == 2) then

    write(*,*) 'Huckel guess...'
    call huckel_guess(nBas,S,Hc,X,c)

  elseif(guess_type == 3) then

    call random_number(c)

    write(*,*) 'Random guess...'
    c(:,:) = 2d0*c(:,:) - 1d0

  else

    print*,'Wrong guess option'
    stop

  end if

end subroutine
