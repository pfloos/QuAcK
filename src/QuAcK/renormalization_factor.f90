subroutine renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,e,Omega,rho,rhox,Z)

! Compute renormalization factor for GW

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: SOSEX
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas),Omega(nS),rho(nBas,nBas,nS),rhox(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,x,jb
  double precision              :: eps
  double precision,allocatable  :: SigC(:),dSigC(:),d2SigC(:)
  double precision,external     :: Z_dcgw

! Output variables

  double precision,intent(out)  :: Z(nBas)

! Allocate

  allocate(SigC(nBas),dSigC(nBas),d2SigC(nBas))

  SigC(:)   = 0d0
  dSigC(:)  = 0d0
  d2SigC(:) = 0d0

! Occupied part of the correlation self-energy

  do x=nC+1,nBas-nR
    do i=nC+1,nO
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
          eps = e(x) - e(i) + Omega(jb)
!         Z(x) = Z(x) + 2d0*Z_dcgw(eps,rho(x,i,jb))
!         SigC(x)   = SigC(x)   + 2d0*rho(x,i,jb)**2/eps
          dSigC(x)  = dSigC(x)  - 2d0*rho(x,i,jb)**2/eps**2
!         d2SigC(x) = d2SigC(x) + 4d0*rho(x,i,jb)**2/eps**3
        enddo
      enddo
    enddo
  enddo

! Virtual part of the correlation self-energy

  do x=nC+1,nBas-nR
    do a=nO+1,nBas-nR
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
          eps = e(x) - e(a) - Omega(jb)
!         Z(x) = Z(x) + 2d0*Z_dcgw(eps,rho(x,a,jb))
!         SigC(x)   = SigC(x)   + 2d0*rho(x,a,jb)**2/eps
          dSigC(x)  = dSigC(x)  - 2d0*rho(x,a,jb)**2/eps**2
!         d2SigC(x) = d2SigC(x) + 4d0*rho(x,a,jb)**2/eps**3
        enddo
      enddo
    enddo
  enddo

  ! SOSEX correction

  if(SOSEX) then

    ! Occupied part of the correlation self-energy

    do x=nC+1,nBas-nR
      do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = e(x) - e(i) + Omega(jb)
            dSigC(x) = dSigC(x) - (rho(x,i,jb)/eps)*(rhox(x,i,jb)/eps)
          enddo
        enddo
      enddo
    enddo

    ! Virtual part of the correlation self-energy

    do x=nC+1,nBas-nR
      do a=nO+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = e(x) - e(a) - Omega(jb)
            dSigC(x) = dSigC(x) - (rho(x,a,jb)/eps)*(rhox(x,a,jb)/eps)
          enddo
        enddo
      enddo
    enddo

  endif

! Compute renormalization factor from derivative of SigC
 
  Z(:) = 1d0/(1d0-dSigC(:))

! Z(:) = 1d0 - dSigC(:) + sqrt( (1d0 - dSigC(:))**2 - 2d0*SigC(:)*d2SigC(:) )
! Z(:) = Z(:)/(SigC(:)*d2SigC(:))

end subroutine renormalization_factor
