subroutine R_irred_Parquet_self_energy(nOrb,nC,nO,nV,nR,e,EcGM,SigC,Z)

! Compute correlation part of the self-energy with only irreducible vertices contribution
  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  double precision,intent(in)   :: e(nOrb)
! Local variables
  integer                       :: p,i,j,a,b
  double precision              :: D2p1h,D2h1p
! Output variables
  double precision,intent(out)  :: EcGM
  double precision,intent(out)  :: SigC(nOrb)
  double precision,intent(out)  :: Z(nOrb)

!----------------------------!
! Static Parquet self-energy !
!----------------------------!
  SigC(:) = 0d0
  ! 2h1p part of the correlation self-energy
  do p=nC+1,nOrb-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nOrb-nR

              D2h1p = e(p) + e(a) - e(i) - e(j)
              SigC(p) = SigC(p) !+ 2d0*rho(p,i,m)**2*(1d0-exp(-2d0*s*Dpim*Dpim))/Dpim
              
           end do
        end do
     end do
  end do
  ! 2p1h part of the correlation self-energy
  do p=nC+1,nOrb-nR
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
 
              D2p1h = e(p) + e(i) - e(a) - e(b)
              SigC(p) = SigC(p) !+ 2d0*rho(p,a,m)**2*(1d0-exp(-2d0*s*Dpam*Dpam))/Dpam
              
            end do
         end do
      end do
   end do 
!------------------------!
! Renormalization factor !
!------------------------!
  Z(:)  = 0d0
  ! 2h1p part of the renormlization factor
  do p=nC+1,nOrb-nR
     do i=nC+1,nO
        do j=nC+1,nO
           do a=nO+1,nOrb-nR

              D2h1p = e(p) + e(a) - e(i) - e(j)
              Z(p) = Z(p) 
              
           end do
        end do
     end do
  end do
  ! 2p1h part of the renormlization factor
  do p=nC+1,nOrb-nR
     do i=nC+1,nO
        do a=nO+1,nOrb-nR
           do b=nO+1,nOrb-nR
 
              D2p1h = e(p) + e(i) - e(a) - e(b)
              Z(p) = Z(p) 
              
            end do
         end do
      end do
   end do 

  Z(:) = 1d0/(1d0 - Z(:))

!-------------------------------------!
! Galitskii-Migdal correlation energy !
!-------------------------------------!

  EcGM = 0d0
  ! do i=nC+1,nO
  !   do a=nO+1,nOrb-nR
  !     do m=1,nS

  !     end do
  !   end do
  ! end do

end subroutine R_irred_Parquet_self_energy
