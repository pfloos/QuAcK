subroutine MP3(nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,EHF,e)

! Perform third-order Moller-Plesset calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBasin
  integer,intent(in)            :: nCin
  integer,intent(in)            :: nOin
  integer,intent(in)            :: nVin
  integer,intent(in)            :: nRin
  double precision,intent(in)   :: ENuc,EHF
  double precision,intent(in)   :: e(nBasin)
  double precision,intent(in)   :: ERI(nBasin,nBasin,nBasin,nBasin)

! Local variables

  double precision              :: eps,E2,EcMP2
  double precision              :: eps1,eps2,E3a,E3b,E3c
  double precision              :: EcMP3

  integer                       :: nBas
  integer                       :: nC
  integer                       :: nO
  integer                       :: nV
  integer                       :: nR
  integer                       :: i,j,k,l,a,b,c,d

  double precision,allocatable  :: se(:)
  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)

  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: VVOO(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)


! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Moller-Plesset third-order calculation    |'
  write(*,*)'************************************************'
  write(*,*)

! Spatial to spin orbitals

  nBas = 2*nBasin
  nC   = 2*nCin
  nO   = 2*nOin
  nV   = 2*nVin
  nR   = 2*nRin

  allocate(se(nBas),sERI(nBas,nBas,nBas,nBas))

  call spatial_to_spin_MO_energy(nBasin,e,nBas,se)
  call spatial_to_spin_ERI(nBasin,ERI,nBas,sERI)

! Antysymmetrize ERIs

  allocate(dbERI(nBas,nBas,nBas,nBas))

  call antisymmetrize_ERI(2,nBas,sERI,dbERI)

  deallocate(sERI)

! Form energy denominator

  allocate(eO(nO),eV(nV))

  eO(:) = se(1:nO)
  eV(:) = se(nO+1:nBas)

  deallocate(se)

! Create integral batches

  allocate(OOOO(nO,nO,nO,nO),OOVV(nO,nO,nV,nV),OVVO(nO,nV,nV,nO),VVOO(nV,nV,nO,nO),VVVV(nV,nV,nV,nV))

  OOOO(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,   1:nO  ,   1:nO  )
  OOVV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVVO(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,nO+1:nBas,   1:nO  )
  VVOO(:,:,:,:) = dbERI(nO+1:nBas,nO+1:nBas,   1:nO  ,   1:nO  )
  VVVV(:,:,:,:) = dbERI(nO+1:nBas,nO+1:nBas,nO+1:nBas,nO+1:nBas)

  deallocate(dbERI)

! Compute MP2 energy

  E2 = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR

          eps = eO(i) + eO(j) - eV(a) - eV(b)
         
          E2 = E2 + OOVV(i,j,a,b)*OOVV(i,j,a,b)/eps

        enddo
      enddo
    enddo
  enddo

  EcMP2 = 0.25d0*E2

! Compute MP3 energy

  E3a = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do k=nC+1,nO
        do l=nC+1,nO
          do a=1,nV-nR
            do b=1,nV-nR

              eps1 = eO(i) + eO(j) - eV(a) - eV(b) 
              eps2 = eO(k) + eO(l) - eV(a) - eV(b) 
             
              E3a = E3a + OOVV(i,j,a,b)*OOOO(k,l,i,j)*VVOO(a,b,k,l)/(eps1*eps2)

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  E3b = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do a=1,nV-nR
        do b=1,nV-nR
          do c=1,nV-nR
            do d=1,nV-nR

              eps1 = eO(i) + eO(j) - eV(a) - eV(b) 
              eps2 = eO(i) + eO(j) - eV(c) - eV(d) 
             
              E3b = E3b + OOVV(i,j,a,b)*VVVV(a,b,c,d)*VVOO(c,d,i,j)/(eps1*eps2)

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  E3c = 0d0

  do i=nC+1,nO
    do j=nC+1,nO
      do k=nC+1,nO
        do a=1,nV-nR
          do b=1,nV-nR
            do c=1,nV-nR

              eps1 = eO(i) + eO(j) - eV(a) - eV(b) 
              eps2 = eO(i) + eO(k) - eV(a) - eV(c) 
    
              E3c = E3c + OOVV(i,j,a,b)*OVVO(k,b,c,j)*VVOO(a,c,i,k)/(eps1*eps2)

            enddo
          enddo
        enddo
      enddo
    enddo
  enddo

  EcMP3 = 0.125d0*E3a + 0.125d0*E3b + E3c

  write(*,*)
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32)')           ' MP3 calculation       '
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP2 contribution      ',EcMP2
  write(*,'(A32,1X,F16.10)') ' MP3 contribution      ',EcMP3
  write(*,'(A32)')           '-----------------------'
  write(*,'(A32,1X,F16.10)') ' MP3 correlation energy',             EcMP2 + EcMP3
  write(*,'(A32,1X,F16.10)') ' MP3 total       energy',ENuc + EHF + EcMP2 + EcMP3
  write(*,'(A32)')           '-----------------------'
  write(*,*)

end subroutine MP3
