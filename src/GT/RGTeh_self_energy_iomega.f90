subroutine RGTeh_self_energy_iomega(eta,wcoord,nBas,nC,nO,nV,nR,nS,e,Om,rhoL,rhoR,EcGM,Sig)

! Compute correlation part of the self-energy for GTeh and the renormalization factor

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rhoL(nBas,nBas,nS)
  double precision,intent(in)   :: rhoR(nBas,nBas,nS)
  complex*16,intent(in)         :: wcoord

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,r
  integer                       :: m
  double precision              :: chem_pot
  complex*16                    :: num,eps

! Output variables

  double precision,intent(out)  :: EcGM
  complex*16,intent(out)        :: Sig(nBas,nBas)

! Initialize 

  chem_pot = 0.5d0*(e(nO)+e(nO+1))
  Sig(:,:) = czero
  EcGM     = 0d0

!----------------!
! GW self-energy !
!----------------!

  ! Occupied part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(Sig,wcoord,chem_pot,rhoL,rhoR,nS,nC,nO,nBas,nR,e,Om) &
!$OMP PRIVATE(m,i,q,p,num,eps) &
!$OMP DEFAULT(NONE)
!$OMP DO
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        do m=1,nS
           do i=nC+1,nO

              eps = wcoord - e(i) + chem_pot + Om(m)
              num = rhoL(i,p,m)*rhoR(i,q,m)
              Sig(p,q) = Sig(p,q) + num/eps

           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

! Virtual part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(Sig,wcoord,chem_pot,rhoL,rhoR,nS,nC,nO,nBas,nR,e,Om) &
!$OMP PRIVATE(m,a,q,p,num,eps) &
!$OMP DEFAULT(NONE)
!$OMP DO  
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        do m=1,nS
           do a=nO+1,nBas-nR

              eps = wcoord - e(a) + chem_pot - Om(m)
              num = rhoL(p,a,m)*rhoR(q,a,m)
              Sig(p,q) = Sig(p,q) + num/eps

           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

 ! Galitskii-Migdal correlation energy
  
  EcGM = 0d0
  do m=1,nS
    do a=nO+1,nBas-nR
      do i=nC+1,nO
        eps = e(a) - e(i) + Om(m)
        num = rhoL(i,a,m)*rhoR(i,a,m)
        EcGM = EcGM - real(num*eps/(eps**2 + eta**2))
      end do
    end do
  end do

end subroutine 
