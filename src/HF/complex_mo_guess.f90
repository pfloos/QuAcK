
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

! Local variables  
  integer                       :: iorb
  double precision,allocatable  :: real_c(:,:)

! Output variables

  complex*16,intent(inout)      :: c(nBas,nOrb)

  if(guess_type == 0) then

    write(*,*) 'Reading HF coefficients...'
    c(:,:) = c(:,:)

  elseif(guess_type == 1) then

    write(*,*) 'Core guess...'
    call complex_core_guess(nBas, nOrb, Hc, X, c)
  
  elseif(guess_type == 3) then
    
    allocate(real_c(nBas,nOrb))
    
    call random_number(real_c)

    write(*,*) 'Random guess...'
    c(:,:) = cmplx(2d0*real_c(:,:) - 1d0,0d0,kind=8)
    
    deallocate(real_c)
  
  elseif(guess_type == 2) then

    write(*,*) 'Huckel guess...'
    call complex_huckel_guess(nBas, nOrb,S,Hc, X, c)
  elseif(guess_type == 4) then

    write(*,*) 'Site basis guess for Hubbard...'
    c(:,:)=(0d0,0d0)
    do iorb=1,nOrb
     c(iorb,iorb) = (1d0,0d0)
    enddo
  else
    print*,'Wrong guess option'
    stop

  end if

end subroutine
