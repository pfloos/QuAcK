subroutine initialize_random_generator(iSeed)

! Initialize random number generator

  implicit none

! Input variables

  integer,intent(in)            :: iSeed

! Local variables

  integer,allocatable           :: Seed(:)
  integer                       :: nSeed

  call random_seed(size = nSeed)
  allocate(Seed(nSeed))
  call random_seed(get=Seed)
  if(iSeed /= 0) then
    Seed = 0
    Seed(1) = iSeed
  endif
  call random_seed(put=Seed)

end subroutine initialize_random_generator
