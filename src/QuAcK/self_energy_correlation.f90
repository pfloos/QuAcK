subroutine self_energy_correlation(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,e,Omega,rho,rhox,EcGM,SigC)

! Compute correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: Omega(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  double precision,intent(in)   :: rhox(nBas,nBas,nS)

! Local variables

  integer                       :: i,j,a,b,p,x,y,jb
  double precision              :: eps

! Output variables

  double precision,intent(out)  :: SigC(nBas,nBas)
  double precision,intent(out)  :: EcGM

! Initialize 

  SigC = 0d0

! COHSEX static approximation

  if(COHSEX) then

   ! COHSEX: SEX of the COHSEX correlation self-energy

    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        do i=nC+1,nO
          do jb=1,nS
            SigC(x,y) = SigC(x,y) + 4d0*rho(x,i,jb)*rho(y,i,jb)/Omega(jb)
          enddo
        enddo
      enddo
    enddo
 
    ! COHSEX: COH part of the COHSEX correlation self-energy
 
    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        do p=nC+1,nBas-nR
          do jb=1,nS
            SigC(x,y) = SigC(x,y) - 2d0*rho(x,p,jb)*rho(y,p,jb)/Omega(jb)
          enddo
        enddo
      enddo
    enddo

    EcGM = 0d0
    do i=nC+1,nO
      EcGM = EcGM + 0.5d0*SigC(i,i)
    enddo

  else

  ! Occupied part of the correlation self-energy

    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        do i=nC+1,nO
          do jb=1,nS
            eps = e(x) - e(i) + Omega(jb)
            SigC(x,y) = SigC(x,y) + 2d0*rho(x,i,jb)*rho(y,i,jb)*eps/(eps**2 + eta**2)
          enddo
        enddo
      enddo
    enddo
 
    ! Virtual part of the correlation self-energy
 
    do x=nC+1,nBas-nR
      do y=nC+1,nBas-nR
        do a=nO+1,nBas-nR
          do jb=1,nS
            eps = e(x) - e(a) - Omega(jb)
            SigC(x,y) = SigC(x,y) + 2d0*rho(x,a,jb)*rho(y,a,jb)*eps/(eps**2 + eta**2)
          enddo
        enddo
      enddo
    enddo

    if(SOSEX) then

      ! SOSEX: occupied part of the correlation self-energy

       do x=nC+1,nBas-nR
        do y=nC+1,nBas-nR
          do i=nC+1,nO
            do jb=1,nS
              eps = e(x) - e(i) + Omega(jb) 
              SigC(x,y) = SigC(x,y) - rho(x,i,jb)*rhox(y,i,jb)*eps/(eps**2 + eta**2)
            enddo
          enddo
        enddo
      enddo
  
      ! SOSEX: virtual part of the correlation self-energy
  
      do x=nC+1,nBas-nR
        do y=nC+1,nBas-nR
          do a=nO+1,nBas-nR
            do jb=1,nS
              eps = e(x) - e(a) - Omega(jb) 
              SigC(x,y) = SigC(x,y) - rho(x,a,jb)*rhox(y,a,jb)*eps/(eps**2 + eta**2)
            enddo
          enddo
        enddo
      enddo

    endif

  endif

end subroutine self_energy_correlation
