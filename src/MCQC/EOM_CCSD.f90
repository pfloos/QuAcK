subroutine EOM_CCSD(ispin,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,e,ERI)

! Compute EOM-CCSD excitation energies: see Stanton & Bartlett JCP 98 7029 (1993)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: nH,nP,nHH,nPP,nSCF,n_diis
  double precision              :: Conv
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: B_ADC(:,:),X_ADC(:,:),e_ADC(:),SigInf(:,:),G_ADC(:,:)
  double precision,allocatable  :: db_ERI(:,:,:,:),eOld(:),error_diis(:,:),e_diis(:,:)

  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: p,q,r,s
  integer                       :: nADC,iADC,jADC


! Hello world

  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'|      EOM-CCSD calculation       |'
  write(*,*)'***********************************'
  write(*,*)

! Number of holes 

  nH  = nO
  nHH = nH*nH

! Number of particles

  nP  = nV
  nPP = nP*nP

  write(*,*) 'Total    states: ',nH + nP
  write(*,*) 'Hole     states: ',nH
  write(*,*) 'Particle states: ',nP

! Size of EOM-CCSD matrices

  nEOM = nH + nP + nH*nPP + nHH*nP + nHH*nPP
  write(*,'(1X,A25,I3,A6,I6)') 'Size of EOM-CCSD matrix: ',nEOM,' x ',nEOM

! Memory allocation

  allocate()

! Construct EOM-CCSD matrix

  H(:,:) = 0d0

  iEOM = 1
  jEOM = 1

  H(iEOM,jEOM) = ECCSD

  do p=1,nO

    jADC = jADC + 1
    B_ADC(jADC,jADC) = e(p)

  enddo

end subroutine EOM_CCSD
