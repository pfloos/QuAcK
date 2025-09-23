subroutine RGW_self_energy_iomega(wcoord,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,Sig)

! Compute correlation part of the self-energy and the renormalization factor 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)
  complex*16,intent(in)         :: wcoord

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,m
  double precision              :: chem_pot
  complex*16                    :: num,eps

! Output variables

  complex*16,intent(out)        :: Sig(nOrb,nOrb)

!----------------!
! GW self-energy !
!----------------!

  chem_pot = 0.5d0*(e(nO)+e(nO+1))
  Sig(:,:) = czero

! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,rho,wcoord,chem_pot,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do i=nC+1,nO
 
          eps = wcoord - e(i) + chem_pot + Om(m)
          num = 2d0*rho(p,i,m)*rho(q,i,m)
          Sig(p,q) = Sig(p,q) + num/eps
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,rho,wcoord,chem_pot,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,a,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do a=nO+1,nOrb-nR
 
          eps = wcoord - e(a) + chem_pot - Om(m)
          num = 2d0*rho(p,a,m)*rho(q,a,m)
          Sig(p,q) = Sig(p,q) + num/eps
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

end subroutine 
