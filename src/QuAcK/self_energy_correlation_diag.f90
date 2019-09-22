subroutine self_energy_correlation_diag(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,rhox,EcGM,SigC)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
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
  double precision,intent(in)   :: rhox(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,p,x,jb
  double precision              :: eps
  double precision,external     :: SigC_dcgw

! Output variables

  double precision,intent(out)  :: SigC(nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC = 0d0

! COHSEX static approximation

  if(COHSEX) then

    ! COHSEX: SEX part of the COHSEX correlation self-energy

    do x=nC+1,nBas-nR
      do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            SigC(x) = SigC(x) + 4d0*rho(x,i,jb)**2/Omega(jb)
          enddo
        enddo
      enddo
    enddo
 
    ! COHSEX: COH part of the COHSEX correlation self-energy
 
    do x=nC+1,nBas-nR
      do p=nC+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            SigC(x) = SigC(x) - 2d0*rho(x,p,jb)**2/Omega(jb)
          enddo
        enddo
      enddo
    enddo

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      EcGM = EcGM + 0.5d0*SigC(i)
    enddo
 
  else
 
    ! Occupied part of the correlation self-energy
 
    do x=nC+1,nBas-nR
      do i=nC+1,nO
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = e(x) - e(i) + Omega(jb)
!           SigC(x) = SigC(x) + 4d0*rho(x,i,jb)**2/(eps + eps*sqrt(1d0 + rho(x,i,jb)**2/eps**2))
            SigC(x) = SigC(x) + 2d0*rho(x,i,jb)**2*eps/(eps**2 + eta**2)
!           SigC(x) = SigC(x) + 2d0*SigC_dcgw(eps,rho(x,i,jb))
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
!           SigC(x) = SigC(x) + 4d0*rho(x,a,jb)**2/(eps + eps*sqrt(1d0 + 4d0*rho(x,a,jb)**2/eps**2))
            SigC(x) = SigC(x) + 2d0*rho(x,a,jb)**2*eps/(eps**2 + eta**2)
!           SigC(x) = SigC(x) + 2d0*SigC_dcgw(eps,rho(x,a,jb))
          enddo
        enddo
      enddo
    enddo

    ! GM correlation energy

    EcGM = 0d0
    do i=nC+1,nO
      do a=nO+1,nBas-nR
        jb = 0
        do j=nC+1,nO
          do b=nO+1,nBas-nR
            jb = jb + 1
            eps = e(a) - e(i) + Omega(jb)
            EcGM = EcGM - 2d0*rho(a,i,jb)*rho(a,i,jb)*eps/(eps**2 + eta**2)
          enddo
        enddo
      enddo
    enddo

    if(SOSEX) then

      ! SOSEX: occupied part of the correlation self-energy

      do x=nC+1,nBas-nR
        do i=nC+1,nO
          jb = 0
          do j=nC+1,nO
            do b=nO+1,nBas-nR
              jb = jb + 1
              eps = e(x) - e(i) + Omega(jb)
              SigC(x) = SigC(x) - rho(x,i,jb)*rhox(x,i,jb)*eps/(eps**2 + eta**2)
            enddo
          enddo
        enddo
      enddo

      ! SOSEX: virtual part of the correlation self-energy

      do x=nC+1,nBas-nR
        do a=nO+1,nBas-nR
          jb = 0
          do j=nC+1,nO
            do b=nO+1,nBas-nR
              jb = jb + 1
              eps = e(x) - e(a) - Omega(jb)
              SigC(x) = SigC(x) - rho(x,a,jb)*rhox(x,a,jb)*eps/(eps**2 + eta**2)
            enddo
          enddo
        enddo
      enddo

      ! GM correlation energy

      do i=nC+1,nO
        do a=nO+1,nBas-nR
          jb = 0
          do j=nC+1,nO
            do b=nO+1,nBas-nR
              jb = jb + 1
              eps = e(a) - e(i) + Omega(jb)
              EcGM = EcGM + rho(a,i,jb)*rhox(a,i,jb)*eps/(eps**2 + eta**2)
            enddo
          enddo
        enddo
      enddo

    endif

  endif

end subroutine self_energy_correlation_diag
