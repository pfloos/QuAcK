subroutine print_GHF(nBas,nBas2,nO,e,C,P,ENuc,ET,EV,EJ,EK,EHF,dipole)

! Print one-electron energies and other stuff for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nBas2
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: e(nBas2)
  double precision,intent(in)        :: C(nBas2,nBas2)
  double precision,intent(in)        :: P(nBas2,nBas2)
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET
  double precision,intent(in)        :: EV
  double precision,intent(in)        :: EJ
  double precision,intent(in)        :: EK
  double precision,intent(in)        :: EHF
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: ixyz
  integer                            :: mu,nu
  integer                            :: HOMO
  integer                            :: LUMO
  double precision                   :: Gap
  double precision                   :: Sx2,Sy2,Sz2,S2

  double precision,allocatable       :: Ca(:,:)
  double precision,allocatable       :: Cb(:,:)
  double precision,allocatable       :: Paa(:,:)
  double precision,allocatable       :: Pab(:,:)
  double precision,allocatable       :: Pba(:,:)
  double precision,allocatable       :: Pbb(:,:)

  double precision,external     :: trace_matrix

! HOMO and LUMO

  HOMO = nO
  LUMO = HOMO + 1
  Gap = e(LUMO)-e(HOMO)

! Density matrices

  allocate(Paa(nBas2,nBas2),Pab(nBas2,nBas2),Pba(nBas2,nBas2),Pbb(nBas2,nBas2))

  Paa(:,:) = P(     1:nBas ,     1:nBas )
  Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
  Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
  Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

  allocate(Ca(nBas,nBas2),Cb(nBas,nBas2))

  Ca(:,:) = C(     1:nBas ,1:nBas2)
  Cb(:,:) = C(nBas+1:nBas2,1:nBas2)

! Compute expectation values of S^2

  Sx2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab+Pba)**2
  do mu=1,nBas
    do nu=1,nBas
        Sx2 = Sx2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) + Pab(mu,nu)*Pab(nu,mu))
    end do
  end do

  Sy2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) - 0.25d0*trace_matrix(nBas,Pab+Pba)**2
  do mu=1,nBas
    do nu=1,nBas
        Sy2 = Sy2 - 0.5d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
    end do
  end do

  Sz2 = 0.25d0*trace_matrix(nBas,Paa+Pbb) + 0.25d0*trace_matrix(nBas,Pab-Pba)**2
  do mu=1,nBas
    do nu=1,nBas
        Sz2 = Sz2 - 0.25d0*(Paa(mu,nu)*Pbb(nu,mu) - Pab(mu,nu)*Pab(nu,mu))
        Sz2 = Sz2 + 0.25d0*(Pab(mu,nu)*Pba(nu,mu) - Pba(mu,nu)*Pab(nu,mu))
    end do
  end do

  print*,Sx2,Sy2,Sz2
  
  S2 = Sx2 + Sy2 + Sz2

! Dump results

  write(*,*)
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32)')           ' Summary              '
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' One-electron energy: ',ET + EV,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Kinetic      energy: ',ET,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Potential    energy: ',EV,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' Two-electron energy: ',EJ + EK,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Hartree      energy: ',EJ,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Exchange     energy: ',EK,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.10,A3)') ' Electronic   energy: ',EHF,' au'
  write(*,'(A32,1X,F16.10,A3)') ' Nuclear   repulsion: ',ENuc,' au'
  write(*,'(A32,1X,F16.10,A3)') ' GHF          energy: ',EHF + ENuc,' au'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.6,A3)')  ' GHF HOMO     energy: ',e(HOMO)*HaToeV,' eV'
  write(*,'(A32,1X,F16.6,A3)')  ' GHF LUMO     energy: ',e(LUMO)*HaToeV,' eV'
  write(*,'(A32,1X,F16.6,A3)')  ' GHF HOMO-LUMO gap  : ',Gap*HaToeV,' eV'
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A32,1X,F16.6)')     ' <S**2>             :',S2
  write(*,'(A50)')           '-----------------------------------------'
  write(*,'(A35)')           ' Dipole moment (Debye)    '
  write(*,'(10X,4A10)')      'X','Y','Z','Tot.'
  write(*,'(10X,4F10.6)')    (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
  write(*,'(A50)')           '-----------------------------------------'
  write(*,*)

! Print results

  write(*,'(A50)')  '---------------------------------------'
  write(*,'(A32)') 'MO coefficients'
  write(*,'(A50)')  '---------------------------------------'
  call matout(nBas2,nBas2,c)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A32)') 'MO energies'
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas2,1,e)
  write(*,*)

end subroutine 
