
! ---

subroutine mo_guess(nBas_AOs, nBas_MOs, guess_type, S, Hc, X, c)

!  Guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs)

! Output variables

  double precision,intent(inout)  :: c(nBas_AOs,nBas_MOs)

  if(guess_type == 0) then

    write(*,*) 'Reading HF coefficients...'
    c(:,:) = c(:,:)

  elseif(guess_type == 1) then

    write(*,*) 'Core guess...'
    call core_guess(nBas_AOs, nBas_MOs, Hc, X, c)

  elseif(guess_type == 2) then

    write(*,*) 'Huckel guess...'
    call huckel_guess(nBas_AOs, nBas_MOs, S, Hc, X, c)

  elseif(guess_type == 3) then

    call random_number(c)

    write(*,*) 'Random guess...'
    c(:,:) = 2d0*c(:,:) - 1d0

  else

    print*,'Wrong guess option'
    stop

  end if

end subroutine
