
! ---

subroutine mo_guess(nBas, nOrb, guess_type, S, Hc, X, c)

!  Guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)

! Local variables

  integer                       :: iorb

! Output variables

  double precision,intent(inout)  :: c(nBas,nOrb)

  if(guess_type == 0) then

    write(*,*) 'Reading HF coefficients...'
    c(:,:) = c(:,:)

  elseif(guess_type == 1) then

    write(*,*) 'Core guess...'
    call core_guess(nBas, nOrb, Hc, X, c)

  elseif(guess_type == 2) then

    write(*,*) 'Huckel guess...'
    call huckel_guess(nBas, nOrb, S, Hc, X, c)

  elseif(guess_type == 3) then

    call random_number(c)

    write(*,*) 'Random guess...'
    c(:,:) = 2d0*c(:,:) - 1d0

  elseif(guess_type == 4) then

    write(*,*) 'Site basis guess for Hubbard...'
    c(:,:)=0d0
    do iorb=1,nOrb
     c(iorb,iorb) = 1d0
    enddo

  else

    print*,'Wrong guess option'
    stop

  end if

end subroutine
