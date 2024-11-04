subroutine GGW_self_energy(eta,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,Sig,Z)

! Compute correlation part of the self-energy and the renormalization factor 

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: e(nOrb)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nOrb,nOrb,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,m
  double precision              :: num,eps

! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: Sig(nOrb,nOrb)
  double precision,intent(out)  :: Z(nOrb)

!----------------!
! GW self-energy !
!----------------!

  Sig(:,:) = 0d0

! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do i=nC+1,nO
 
          eps = e(p) - e(i) + Om(m)
          num = rho(p,i,m)*rho(q,i,m)
          Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,a,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do q=nC+1,nOrb-nR
    do p=nC+1,nOrb-nR
      do m=1,nS
        do a=nO+1,nOrb-nR
 
          eps = e(p) - e(a) - Om(m)
          num = rho(p,a,m)*rho(q,a,m)
          Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
 
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

!------------------------!
! Renormalization factor !
!------------------------!

  Z(:)     = 0d0

! Occupied part of the renormalization factor

!$OMP PARALLEL &
!$OMP SHARED(Sig,Z,rho,eta,nS,nC,nO,nOrb,nR,e,Om) &
!$OMP PRIVATE(m,i,q,p,eps,num) &
!$OMP DEFAULT(NONE)
!$OMP DO
  do p=nC+1,nOrb-nR
     do m=1,nS
        do i=nC+1,nO
 
           eps = e(p) - e(i) + Om(m)
           num = rho(p,i,m)*rho(q,i,m)
           Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
 
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the renormalization factor

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nS,nC,nO,nOrb,nR,e,Om) &
  !$OMP PRIVATE(m,a,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO  
  do p=nC+1,nOrb-nR
    do m=1,nS
      do a=nO+1,nOrb-nR
 
        eps = e(p) - e(a) - Om(m)
        num = rho(p,a,m)*rho(q,a,m)
        Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
 
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  Z(:) = 1d0/(1d0 - Z(:))

!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  do m=1,nS
    do a=nO+1,nOrb-nR
      do i=nC+1,nO

        eps = e(a) - e(i) + Om(m)
        num = rho(a,i,m)*rho(a,i,m)
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

end subroutine 
