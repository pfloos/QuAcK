subroutine GW_self_energy(COHSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,EcGM,SigC)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
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
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nBas,nBas)

! Initialize 

  SigC(:,:) = 0d0

!-----------------------------!
! COHSEX static approximation !
!-----------------------------!

  if(COHSEX) then

   ! COHSEX: SEX of the COHSEX correlation self-energy

    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do i=nC+1,nO
          do jb=1,nS
            SigC(p,q) = SigC(p,q) + 4d0*rho(p,i,jb)*rho(q,i,jb)/Omega(jb)
          end do
        end do
      end do
    end do
 
    ! COHSEX: COH part of the COHSEX correlation self-energy
 
    do p=nC+1,nBas-nR
      do q=nC+1,nBas-nR
        do r=nC+1,nBas-nR
          do jb=1,nS
            SigC(p,q) = SigC(p,q) - 2d0*rho(p,r,jb)*rho(q,r,jb)/Omega(jb)
          end do
        end do
      end do
    end do

    EcGM = 0d0
    do i=nC+1,nO
      EcGM = EcGM + 0.5d0*SigC(i,i)
    end do

  else

!----------------!
! GW self-energy !
!----------------!

  ! Occupied part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(SigC,rho,eta,nS,nC,nO,nBas,nR,e,Omega) &
!$OMP PRIVATE(jb,i,q,p,eps) &
!$OMP DEFAULT(NONE)
!$OMP DO
do q=nC+1,nBas-nR
   do p=nC+1,nBas-nR
      do jb=1,nS
         do i=nC+1,nO
            eps = e(p) - e(i) + Omega(jb)
            SigC(p,q) = SigC(p,q) + 2d0*rho(p,i,jb)*rho(q,i,jb)*eps/(eps**2 + eta**2)
         end do
      end do
   end do
end do
!$OMP END DO
!$OMP END PARALLEL

    ! Virtual part of the correlation self-energy

!$OMP PARALLEL &
!$OMP SHARED(SigC,rho,eta,nS,nC,nO,nBas,nR,e,Omega) &
!$OMP PRIVATE(jb,a,q,p,eps) &
!$OMP DEFAULT(NONE)
!$OMP DO  
do q=nC+1,nBas-nR
   do p=nC+1,nBas-nR
      do jb=1,nS
         do a=nO+1,nBas-nR
            eps = e(p) - e(a) - Omega(jb)
            SigC(p,q) = SigC(p,q) + 2d0*rho(p,a,jb)*rho(q,a,jb)*eps/(eps**2 + eta**2)
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
          EcGM = EcGM - 4d0*rho(a,i,jb)*rho(a,i,jb)*eps/(eps**2 + eta**2)
        end do
      end do
    end do

  end if

end subroutine 
