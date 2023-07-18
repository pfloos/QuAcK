subroutine GW_self_energy(eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,Sig,Z)

! Compute correlation part of the self-energy and the renormalization factor 

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
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b
  integer                       :: p,q,r
  integer                       :: jb
  double precision              :: num,eps

! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: Sig(nBas,nBas)
  double precision,intent(out)  :: Z(nBas)

! Initialize 

  Sig(:,:) = 0d0
  Z(:)     = 0d0

!----------------!
! GW self-energy !
!----------------!

! Occupied part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(Sig,rho,eta,nS,nC,nO,nBas,nR,e,Omega) &
!$OMP PRIVATE(jb,i,q,p,eps,num) &
!$OMP DEFAULT(NONE)
!$OMP DO
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        do jb=1,nS
           do i=nC+1,nO
 
              eps = e(p) - e(i) + Omega(jb)
              num = 2d0*rho(p,i,jb)*rho(q,i,jb)
              Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
              if(p == q) Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
 
           end do
        end do
     end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

! Virtual part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(Sig,rho,eta,nS,nC,nO,nBas,nR,e,Omega) &
!$OMP PRIVATE(jb,a,q,p,eps,num) &
!$OMP DEFAULT(NONE)
!$OMP DO  
  do q=nC+1,nBas-nR
     do p=nC+1,nBas-nR
        do jb=1,nS
           do a=nO+1,nBas-nR
 
              eps = e(p) - e(a) - Omega(jb)
              num = 2d0*rho(p,a,jb)*rho(q,a,jb)
              Sig(p,q) = Sig(p,q) + num*eps/(eps**2 + eta**2)
              if(p == q) Z(p) = Z(p) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
 
           end do
        end do
     end do
  end do
!$OMP END DO
!$OMP END PARALLEL

! Galitskii-Migdal correlation energy

  EcGM = 0d0
  do jb=1,nS
    do a=nO+1,nBas-nR
      do i=nC+1,nO

        eps = e(a) - e(i) + Omega(jb)
        num = 4d0*rho(a,i,jb)*rho(a,i,jb)
        EcGM = EcGM - num*eps/(eps**2 + eta**2)

      end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
