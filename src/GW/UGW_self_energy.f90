subroutine UGW_self_energy(eta,nBas,nC,nO,nV,nR,nSt,e,Om,rho,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: Om(nSt)
  double precision,intent(in)   :: rho(nBas,nBas,nSt,nspin)

! Local variables

  integer                       :: i,a,p,q,m
  double precision              :: num,eps
  double precision              :: EcGM_alpha, EcGM_beta

! Output variables

  double precision,intent(out)  :: Sig(nBas,nBas,nspin)
  double precision,intent(out)  :: Z(nBas,nspin)
  double precision              :: EcGM(nspin)

! Initialize 

  Sig(:,:,:) = 0d0
  Z(:,:)     = 0d0
  EcGM(:)    = 0d0

!--------------!
! Spin-up part !
!--------------!

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC(1)+1,nBas-nR(1)
    do p=nC(1)+1,nBas-nR(1)
      do m=1,nSt
        do i=nC(1)+1,nO(1)
          eps = e(p,1) - e(i,1) + Om(m)
          num = rho(p,i,m,1)*rho(q,i,m,1)
          Sig(p,q,1) = Sig(p,q,1) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,a,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC(1)+1,nBas-nR(1)
    do p=nC(1)+1,nBas-nR(1)
      do m=1,nSt
        do a=nO(1)+1,nBas-nR(1)
          eps = e(p,1) - e(a,1) - Om(m)
          num = rho(p,a,m,1)*rho(q,a,m,1)
          Sig(p,q,1) = Sig(p,q,1) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p,1) = Z(p,1) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! GM correlation energy
  
  EcGM_alpha = 0d0
  !$OMP PARALLEL &
  !$OMP SHARED(rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,i,a,eps,num) &
  !$OMP DEFAULT(NONE) &
  !$OMP REDUCTION(-:EcGM_alpha)
  !$OMP DO
  do m=1,nSt
    do a=nO(1)+1,nBas-nR(1)
      do i=nC(1)+1,nO(1)
        eps = e(a,1) - e(i,1) + Om(m)
        num = rho(a,i,m,1)**2
        EcGM_alpha = EcGM_alpha - num*eps/(eps**2 + eta**2)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  EcGM(1) = EcGM_alpha

!----------------!
! Spin-down part !
!----------------!

  ! Occupied part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,i,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC(2)+1,nBas-nR(2)
    do p=nC(2)+1,nBas-nR(2)
      do m=1,nSt
        do i=nC(2)+1,nO(2)
          eps = e(p,2) - e(i,2) + Om(m)
          num = rho(p,i,m,2)*rho(q,i,m,2)
          Sig(p,q,2) = Sig(p,q,2) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! Virtual part of the correlation self-energy

  !$OMP PARALLEL &
  !$OMP SHARED(Sig,Z,rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,a,q,p,eps,num) &
  !$OMP DEFAULT(NONE)
  !$OMP DO
  do q=nC(2)+1,nBas-nR(2)
    do p=nC(2)+1,nBas-nR(2)
      do m=1,nSt
        do a=nO(2)+1,nBas-nR(2)
          eps = e(p,2) - e(a,2) - Om(m)
          num = rho(p,a,m,2)*rho(q,a,m,2)
          Sig(p,q,2) = Sig(p,q,2) + num*eps/(eps**2 + eta**2)
          if(p == q) Z(p,2) = Z(p,2) - num*(eps**2 - eta**2)/(eps**2 + eta**2)**2
        end do
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL

  ! GM correlation energy
  
  EcGM_beta = 0d0
  !$OMP PARALLEL &
  !$OMP SHARED(rho,eta,nSt,nC,nO,nBas,nR,e,Om) &
  !$OMP PRIVATE(m,i,a,eps,num) &
  !$OMP DEFAULT(NONE) &
  !$OMP REDUCTION(-:EcGM_beta)
  !$OMP DO
  do m=1,nSt
    do a=nO(2)+1,nBas-nR(2)
      do i=nC(2)+1,nO(2)
        eps = e(a,2) - e(i,2) + Om(m)
        num = rho(a,i,m,2)**2
        EcGM_beta = EcGM_beta - num*eps/(eps**2 + eta**2)
      end do
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL
  EcGM(2) = EcGM_beta

! Compute renormalization factor from derivative 

  Z(:,:) = 1d0/(1d0 - Z(:,:))

end subroutine 
