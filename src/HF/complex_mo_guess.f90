
! ---

subroutine complex_mo_guess(nBas, nOrb, guess_type, S, Hc, X, c)

!  Guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: S(nBas,nBas)
  complex*16,intent(in)         :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)

! Output variables

  complex*16,intent(inout)      :: c(nBas,nOrb)

  if(guess_type == 0) then

    write(*,*) 'Reading HF coefficients...'
    c(:,:) = c(:,:)

  elseif(guess_type == 1) then

    write(*,*) 'Core guess...'
    call complex_core_guess(nBas, nOrb, Hc, X, c)
  
  elseif(guess_type == 2) then

    write(*,*) 'Huckel guess...'
    call complex_huckel_guess(nBas, nOrb,S,Hc, X, c)
  else
    print*,'Wrong guess option'
    stop

  end if

end subroutine
