subroutine CID(dotest,singlet,triplet,nBasin,nCin,nOin,nVin,nRin,ERIin,eHFin,E0)

! Perform configuration interaction with doubles

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBasin
  integer,intent(in)            :: nCin
  integer,intent(in)            :: nOin
  integer,intent(in)            :: nVin
  integer,intent(in)            :: nRin
  double precision,intent(in)   :: eHFin(nBasin)
  double precision,intent(in)   :: ERIin(nBasin,nBasin,nBasin,nBasin)
  double precision,intent(in)   :: E0

! Local variables

  integer                       :: nBas
  integer                       :: nC
  integer                       :: nO
  integer                       :: nV
  integer                       :: nR

  double precision,allocatable  :: eHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: ERI(:,:,:,:)

  logical                       :: dump_trans = .false.
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,kc,iajb,kcld
  integer                       :: ishift,jshift
  integer                       :: ispin
  integer                       :: nD
  integer                       :: nH
  integer                       :: maxH
  double precision,external     :: Kronecker_delta
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: ECID(:)

  double precision              :: tmp

! Hello world

  write(*,*)
  write(*,*)'******************************************************'
  write(*,*)'| Configuration Interaction with Singles and Doubles |'
  write(*,*)'******************************************************'
  write(*,*)

! Spatial to spin orbitals

  nBas = 2*nBasin
  nC   = 2*nCin
  nO   = 2*nOin
  nV   = 2*nVin
  nR   = 2*nRin

  allocate(eHF(nBas),sERI(nBas,nBas,nBas,nBas))
                      
  call spatial_to_spin_MO_energy(nBasin,eHFin,nBas,eHF)
  call spatial_to_spin_ERI(nBasin,ERIin,nBas,sERI)

! Antysymmetrize ERIs

  allocate(ERI(nBas,nBas,nBas,nBas))

  call antisymmetrize_ERI(2,nBas,sERI,ERI)

  deallocate(sERI)

! Compute CID matrix

  nD  = (nO - nC)*(nO - nC - 1)/2*(nV - nR)*(nV - nR - 1)/2
  nH  = 1 + nD

  write(*,*) 'nD = ',nD 
  write(*,*) 'nH = ',nH
  write(*,*)

  maxH = min(nH,21)

  ! Memory allocation

  allocate(H(nH,nH),ECID(nH))
 
  ! 00 block

  ishift = 0 
  jshift = 0

  H(ishift+1,jshift+1) = E0
 
  print*,'00 block  done...'

  ! 0D blocks

  ishift = 0
  jshift = 1

  iajb = 0
  do i=nC+1,nO
    do a=1,nV-nR
      do j=i+1,nO
        do b=a+1,nV-nR

          iajb = iajb + 1
          tmp = ERI(i,j,nO+a,nO+b)

          H(ishift+1,jshift+iajb) = tmp
          H(jshift+iajb,ishift+1) = tmp

        end do
      end do
    end do
  end do
  
  print*,'0D blocks done...'

  ! DD block

  ishift = 1 
  jshift = 1

  iajb = 0
  do i=nC+1,nO
    do a=1,nV-nR
      do j=i+1,nO
        do b=a+1,nV-nR

          iajb = iajb + 1

          kcld = 0
          do k=nC+1,nO
            do c=1,nV-nR
              do l=k+1,nO
                do d=c+1,nV-nR
       
                  kcld = kcld + 1
                  tmp = &
                        E0*Kronecker_delta(i,k)*Kronecker_delta(j,l)*Kronecker_delta(a,c)*Kronecker_delta(b,d) & 
                      + eHF(j)*Kronecker_delta(l,j)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(i,k) &
                      - eHF(j)*Kronecker_delta(l,j)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(i,k) &
                      - eHF(j)*Kronecker_delta(k,j)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(i,l) &
                      + eHF(j)*Kronecker_delta(k,j)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(i,l) &
                      - eHF(i)*Kronecker_delta(l,i)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(j,k) &
                      + eHF(i)*Kronecker_delta(l,i)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(j,k) &
                      + eHF(i)*Kronecker_delta(k,i)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(j,l) &
                      - eHF(i)*Kronecker_delta(k,i)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(j,l) &
                      + eHF(a)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(i,l)*Kronecker_delta(j,k) &
                      - eHF(a)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(i,l)*Kronecker_delta(j,k) &
                      - eHF(a)*Kronecker_delta(a,d)*Kronecker_delta(b,c)*Kronecker_delta(i,k)*Kronecker_delta(j,l) &
                      + eHF(a)*Kronecker_delta(a,c)*Kronecker_delta(b,d)*Kronecker_delta(i,k)*Kronecker_delta(j,l) &
                      - eHF(b)*Kronecker_delta(b,d)*Kronecker_delta(a,c)*Kronecker_delta(i,l)*Kronecker_delta(j,k) &
                      + eHF(b)*Kronecker_delta(b,c)*Kronecker_delta(a,d)*Kronecker_delta(i,l)*Kronecker_delta(j,k) &
                      + eHF(b)*Kronecker_delta(b,d)*Kronecker_delta(a,c)*Kronecker_delta(i,k)*Kronecker_delta(j,l) &
                      - eHF(b)*Kronecker_delta(b,c)*Kronecker_delta(a,d)*Kronecker_delta(i,k)*Kronecker_delta(j,l) &  
                      - ERI(k,l,i,j)*Kronecker_delta(a,d)*Kronecker_delta(b,c) &
                      + ERI(k,l,i,j)*Kronecker_delta(a,c)*Kronecker_delta(b,d) &
                      + ERI(nO+a,l,nO+d,j)*Kronecker_delta(b,c)*Kronecker_delta(i,k) & 
                      - ERI(nO+a,l,nO+c,j)*Kronecker_delta(b,d)*Kronecker_delta(i,k) &
                      - ERI(nO+a,k,nO+d,j)*Kronecker_delta(b,c)*Kronecker_delta(i,l) &
                      + ERI(nO+a,k,nO+c,j)*Kronecker_delta(b,d)*Kronecker_delta(i,l) &
                      - ERI(nO+a,l,nO+d,i)*Kronecker_delta(b,c)*Kronecker_delta(j,k) &
                      + ERI(nO+a,l,nO+c,i)*Kronecker_delta(b,d)*Kronecker_delta(j,k) &
                      + ERI(nO+a,k,nO+d,i)*Kronecker_delta(b,c)*Kronecker_delta(j,l) &
                      - ERI(nO+a,k,nO+c,i)*Kronecker_delta(b,d)*Kronecker_delta(j,l) &
                      - ERI(nO+b,l,nO+d,j)*Kronecker_delta(a,c)*Kronecker_delta(i,k) &
                      + ERI(nO+b,l,nO+c,j)*Kronecker_delta(a,d)*Kronecker_delta(i,k) &
                      + ERI(nO+b,k,nO+d,j)*Kronecker_delta(a,c)*Kronecker_delta(i,l) &
                      - ERI(nO+b,k,nO+c,j)*Kronecker_delta(a,d)*Kronecker_delta(i,l) &
                      + ERI(nO+b,l,nO+d,i)*Kronecker_delta(a,c)*Kronecker_delta(j,k) &
                      - ERI(nO+b,l,nO+c,i)*Kronecker_delta(a,d)*Kronecker_delta(j,k) &
                      - ERI(nO+b,k,nO+d,i)*Kronecker_delta(a,c)*Kronecker_delta(j,l) &
                      + ERI(nO+b,k,nO+c,i)*Kronecker_delta(a,d)*Kronecker_delta(j,l) &
                      - ERI(nO+a,nO+b,nO+c,nO+d)*Kronecker_delta(i,l)*Kronecker_delta(j,k) &
                      + ERI(nO+a,nO+b,nO+c,nO+d)*Kronecker_delta(i,k)*Kronecker_delta(j,l)   

                  H(ishift+iajb,jshift+kcld) = tmp
       
                end do
              end do
            end do
          end do

        end do
      end do
    end do
  end do

  print*,'DD block  done...'

  write(*,*)
  write(*,*) 'Diagonalizing CID matrix...'
  write(*,*)

  call diagonalize_matrix(nH,H,ECID)

  print*,'CID energies (au)'
  call matout(maxH,1,ECID)
  write(*,*)

  print*,'CID excitation energies (eV)'
  call matout(maxH-1,1,(ECID(2:maxH)-ECID(1))*HaToeV)
  write(*,*)

  if(dump_trans) then
    print*,'Singlet CID transition vectors'
    call matout(nH,nH,H)
    write(*,*)
  end if

end subroutine 
