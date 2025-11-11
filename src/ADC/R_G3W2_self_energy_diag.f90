subroutine R_G3W2_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,EcGM,Sig,Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor
! for the G3W2 approximation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
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

! Local variables

  integer                       :: p,i,j,k,a,b,c,t,s
  double precision              :: num,dem1,dem2,dem3,dem4

! Output variables

  double precision,intent(out)  :: Sig(nOrb)
  double precision,intent(out)  :: Z(nOrb)
  double precision,intent(out)  :: EcGM

! Initialize 

  Sig(:) = 0d0
  Z(:)   = 0d0
  EcGM   = 0d0

!------------------!
! G3W2 self-energy !
!------------------!
 
! ooo

  do p=nC+1,nOrb-nR
    do i=nC+1,nO
    do j=nC+1,nO
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,i,t)*rho(j,k,t)*rho(p,k,s)*rho(i,j,s)

        dem1 = e(p) - e(i) + Om(t)
        dem2 = e(p) - e(j) + Om(t) + Om(s)
        dem3 = e(p) - e(k) + Om(s)

        Sig(p) = Sig(p) + num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! vvv

  do p=nC+1,nOrb-nR
    do a=nO+1,nOrb-nR
    do b=nO+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,a,t)*rho(b,c,t)*rho(p,c,s)*rho(a,b,s)

        dem1 = e(p) - e(a) - Om(t)
        dem2 = e(p) - e(b) - Om(t) - Om(s)
        dem3 = e(p) - e(c) - Om(s)

        Sig(p) = Sig(p) + num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! voo + oov

  do p=nC+1,nOrb-nR
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,a,t)*rho(j,k,t)*rho(p,k,s)*rho(a,j,s)

        dem1 = e(p) - e(k) + Om(s)
        dem2 = Om(s) + e(a) - e(j)
        dem3 = e(p) - e(a) - Om(t)

        Sig(p) = Sig(p) + 2d0*num/(dem1*dem2*dem3)

        dem1 = e(p) - e(k) + Om(s)
        dem2 = Om(s) + e(a) - e(j)
        dem3 = e(p) - e(j) + Om(s) + Om(t)

        Sig(p) = Sig(p) - 2d0*num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! ovv + vvo

  do p=nC+1,nOrb-nR
    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,i,t)*rho(b,c,t)*rho(p,c,s)*rho(i,b,s)

        dem1 = e(p) - e(c) - Om(s)
        dem2 = Om(s) - e(i) + e(b)
        dem3 = e(p) - e(i) + Om(t)

        Sig(p) = Sig(p) - 2d0*num/(dem1*dem2*dem3)

        dem1 = e(p) - e(c) - Om(s)
        dem2 = Om(s) - e(i) + e(b)
        dem3 = e(p) - e(b) - Om(s) - Om(t)

        Sig(p) = Sig(p) + 2d0*num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! ovo

  do p=nC+1,nOrb-nR
    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,i,t)*rho(b,k,t)*rho(p,k,s)*rho(i,b,s)

        dem1 = e(p) - e(i) - e(k) + e(b)
        dem2 = e(p) - e(b) - Om(t) - Om(s) 
        dem3 = Om(s) + e(b) - e(i) 
        dem4 = Om(t) + e(b) - e(k) 

        Sig(p) = Sig(p) + (2d0*e(b) - e(i) - e(k) + Om(t) + Om(s))*num/(dem1*dem2*dem3*dem4)

        dem1 = e(p) - e(i) - e(k) + e(b)
        dem2 = Om(t) + e(b) - e(k) 
        dem3 = e(p) - e(k) + Om(s) 

        Sig(p) = Sig(p) - 2d0*num/(dem1*dem2*dem3)

        dem1 = e(p) - e(i) - e(k) + e(b)
        dem2 = e(p) - e(i) + Om(t)
        dem3 = e(p) - e(k) + Om(s) 

        Sig(p) = Sig(p) - 1d0*num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! vov

  do p=nC+1,nOrb-nR
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS

        num = 4d0*rho(p,a,t)*rho(j,c,t)*rho(p,c,s)*rho(a,j,s)

        dem1 = e(p) - e(a) - e(c) + e(j)
        dem2 = e(p) - e(j) + Om(t) + Om(s) 
        dem3 = Om(s) - e(j) + e(a) 
        dem4 = Om(t) - e(j) + e(c) 

        Sig(p) = Sig(p) + (2d0*e(j) - e(a) - e(c) - Om(t) - Om(s))*num/(dem1*dem2*dem3*dem4)

        dem1 = e(p) - e(a) - e(c) + e(j)
        dem2 = Om(t) - e(j) + e(c) 
        dem3 = e(p) - e(c) - Om(s) 

        Sig(p) = Sig(p) + 2d0*num/(dem1*dem2*dem3)

        dem1 = e(p) - e(a) - e(c) + e(j)
        dem2 = e(p) - e(a) - Om(t)
        dem3 = e(p) - e(c) - Om(s) 

        Sig(p) = Sig(p) - 1d0*num/(dem1*dem2*dem3)

      end do
      end do
    end do
    end do
    end do
  end do

! Compute renormalization factor from derivative 

  Z(:) = 1d0/(1d0 - Z(:))

end subroutine 
