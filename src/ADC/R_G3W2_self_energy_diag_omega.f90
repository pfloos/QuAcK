subroutine R_G3W2_self_energy_diag_omega(p,w,eta,flow,nBas,nOrb,nC,nO,nV,nR,nS,e,Om,rho,ERI,Sig,Z)

! Compute diagonal of the correlation part of the self-energy and the renormalization factor
! for the G3W2 approximation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: p
  double precision,intent(in)   :: w

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
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(in)   :: flow 

! Local variables

  integer                       :: i,j,k
  integer                       :: a,b,c
  integer                       :: m,t,s
  double precision              :: num,dem,dem1,dem2,dem3,dem4
  double precision              :: reg,reg1,reg2,reg3,reg4

  logical                       :: add_GW
  logical                       :: add_SOX
  logical                       :: add_2SOSEX
  logical                       :: add_G3W2

! Output variables

  double precision,intent(out)  :: Sig
  double precision,intent(out)  :: Z

! Initialize 

  Sig = 0d0
  Z   = 0d0

! Self-energy components

  add_GW     = .true.
  add_SOX    = .true.
  add_2SOSEX = .true.
  add_G3W2   = .true.

!----------------!
! GW self-energy !
!----------------!

  if(add_GW) then

  ! o
 
    do i=nC+1,nO
      do m=1,nS
 
        num = 2d0*rho(p,i,m)*rho(p,i,m)
        dem = w - e(i) + Om(m)
        reg = (1d0 - exp(-2d0*flow*dem*dem))
        Sig = Sig + num*reg/dem
        Z   = Z   - num*reg/(dem*dem)
 
      end do
    end do
 
  ! v
 
    do a=nO+1,nOrb-nR
      do m=1,nS
 
        num = 2d0*rho(p,a,m)*rho(p,a,m)
        dem = w - e(a) - Om(m)
        reg = (1d0 - exp(-2d0*flow*dem*dem))
        Sig = Sig + num*reg/dem
        Z   = Z   - num*reg/(dem*dem)
 
      end do
    end do

  end if

!-----------------!
! SOX self-energy !
!-----------------!

  if(add_SOX) then

  ! ovo

    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do k=nC+1,nO
 
      num = - ERI(p,b,i,k)*ERI(p,b,k,i)
      dem = w - e(i) + e(b) - e(k) 
      reg = (1d0 - exp(-2d0*flow*dem*dem))
      Sig = Sig + num*reg/dem
      Z   = Z   - num*reg/(dem*dem)
 
    end do
    end do
    end do
 
  ! vov

    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do c=nO+1,nOrb-nR
 
      num = - ERI(p,j,a,c)*ERI(p,j,c,a)
      dem = w - e(a) + e(j) - e(c) 
      reg = (1d0 - exp(-2d0*flow*dem*dem))
      Sig = Sig + num*reg/dem
      Z   = Z   - num*reg/(dem*dem)
 
    end do
    end do
    end do

  end if

!-------------------!
! SOSEX self-energy !
!-------------------!

  if(add_2SOSEX) then

  ! ovo

    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do k=nC+1,nO
      do s=1,nS
 
        num = 4d0*ERI(p,b,i,k)*rho(p,k,s)*rho(i,b,s)
       
        dem1 = w - e(i) + e(b) - e(k)
        dem2 = e(b) - e(i) + Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig + num*reg1*reg2/(dem1*dem2)
        Z = Z - num*reg1*reg2/(dem1*dem1*dem2)
       
        dem1 = w - e(i) + e(b) - e(k)
        dem2 = w - e(k) + Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig + num*reg1*reg2/(dem1*dem2)
        Z = Z - num*reg1*reg2/(dem1*dem1*dem2)
       
       
 
      end do
    end do
    end do
    end do
 
  ! vov
 
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do c=nO+1,nOrb-nR
      do s=1,nS
 
        num = 4d0*ERI(p,j,a,c)*rho(p,c,s)*rho(a,j,s)
       
        dem1 = w - e(a) + e(j) - e(c)
        dem2 = e(a) - e(j) + Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig + num*reg1*reg2/(dem1*dem2)
        Z = Z - num*reg1*reg2/(dem1*dem1*dem2)
       
        dem1 = w - e(a) + e(j) - e(c)
        dem2 = w - e(c) - Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig - num*reg1*reg2/(dem1*dem2)
        Z = Z + num*reg1*reg2/(dem1*dem1*dem2)
       
      end do
    end do
    end do
    end do
 
 ! voo
 
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do k=nC+1,nO
      do s=1,nS
 
        num = 4d0*ERI(p,j,a,k)*rho(p,k,s)*rho(a,j,s)
       
        dem1 = e(a) - e(j) + Om(s)
        dem2 = w - e(k) + Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig + num*reg1*reg2/(dem1*dem2)
        Z = Z - num*reg1*reg2/(dem1*dem1*dem2)
       
      end do
    end do
    end do
    end do
 
  ! ovv
 
    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do s=1,nS
 
        num = 4d0*ERI(p,b,i,c)*rho(p,c,s)*rho(i,b,s)
       
        dem1 = e(b) - e(i) + Om(s)
        dem2 = w - e(c) - Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))

        Sig = Sig + num*reg1*reg2/(dem1*dem2)
        Z = Z - num*reg1*reg2/(dem1*dem1*dem2)
       
      end do
    end do
    end do
    end do

  end if

!------------------!
! G3W2 self-energy !
!------------------!
 
  if(add_G3W2) then

  ! ooo

    do i=nC+1,nO
    do j=nC+1,nO
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,i,t)*rho(j,k,t)*rho(p,k,s)*rho(i,j,s)
 
        dem1 = w - e(i) + Om(t)
        dem2 = w - e(j) + Om(t) + Om(s)
        dem3 = w - e(k) + Om(s)
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig + num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z - num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              - num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              - num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)

      end do
      end do
    end do
    end do
    end do
 
  ! vvv
 
    do a=nO+1,nOrb-nR
    do b=nO+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,a,t)*rho(b,c,t)*rho(p,c,s)*rho(a,b,s)
 
        dem1 = w - e(a) - Om(t)
        dem2 = w - e(b) - Om(t) - Om(s)
        dem3 = w - e(c) - Om(s)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig + num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z - num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              - num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              - num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)
 
      end do
      end do
    end do
    end do
    end do
 
  ! voo + oov
 
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,a,t)*rho(j,k,t)*rho(p,k,s)*rho(a,j,s)
 
        dem1 = w - e(k) + Om(s)
        dem2 = Om(s) + e(a) - e(j)
        dem3 = w - e(a) - Om(t)
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z - 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3) 

        dem1 = w - e(k) + Om(s)
        dem2 = Om(s) + e(a) - e(j)
        dem3 = w - e(j) + Om(s) + Om(t)

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z + 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)
 
      end do
      end do
    end do
    end do
    end do
 
  ! ovv + vvo
 
    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,i,t)*rho(b,c,t)*rho(p,c,s)*rho(i,b,s)
 
        dem1 = w - e(c) - Om(s)
        dem2 = Om(s) - e(i) + e(b)
        dem3 = w - e(i) + Om(t)
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z + 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)

        dem1 = w - e(c) - Om(s)
        dem2 = Om(s) - e(i) + e(b)
        dem3 = w - e(b) - Om(s) - Om(t)
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z - 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)

      end do
      end do
    end do
    end do
    end do
 
  ! ovo
 
    do i=nC+1,nO
    do b=nO+1,nOrb-nR
    do k=nC+1,nO
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,i,t)*rho(b,k,t)*rho(p,k,s)*rho(i,b,s)
 
        dem1 = w - e(i) - e(k) + e(b)
        dem2 = w - e(b) - Om(t) - Om(s) 
        dem3 = Om(s) + e(b) - e(i) 
        dem4 = Om(t) + e(b) - e(k) 

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))
        reg4 = (1d0 - exp(-2d0*flow*dem4*dem4))

        Sig = Sig + (2d0*e(b) - e(i) - e(k) + Om(t) + Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem2*dem3*dem4)

        Z = Z - (2d0*e(b) - e(i) - e(k) + Om(t) + Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem1*dem2*dem3*dem4) &
              - (2d0*e(b) - e(i) - e(k) + Om(t) + Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem2*dem2*dem3*dem4)

        dem1 = w - e(i) - e(k) + e(b)
        dem2 = Om(t) + e(b) - e(k) 
        dem3 = w - e(k) + Om(s) 

        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z + 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)
 
        dem1 = w - e(i) - e(k) + e(b)
        dem2 = w - e(i) + Om(t)
        dem3 = w - e(k) + Om(s) 
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig - num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z + num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              + num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              + num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)

 
      end do
      end do
    end do
    end do
    end do
 
  ! vov
 
    do a=nO+1,nOrb-nR
    do j=nC+1,nO
    do c=nO+1,nOrb-nR
      do t=1,nS
      do s=1,nS
 
        num = 4d0*rho(p,a,t)*rho(j,c,t)*rho(p,c,s)*rho(a,j,s)
 
        dem1 = w - e(a) - e(c) + e(j)
        dem2 = w - e(j) + Om(t) + Om(s) 
        dem3 = Om(s) - e(j) + e(a) 
        dem4 = Om(t) - e(j) + e(c) 
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))
        reg4 = (1d0 - exp(-2d0*flow*dem4*dem4))

        Sig = Sig + (2d0*e(j) - e(a) - e(c) - Om(t) - Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem2*dem3*dem4)

        Z = Z - (2d0*e(j) - e(a) - e(c) - Om(t) - Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem1*dem2*dem3*dem4) &
              - (2d0*e(j) - e(a) - e(c) - Om(t) - Om(s))*num*reg1*reg2*reg3*reg4/(dem1*dem2*dem2*dem3*dem4)

        dem1 = w - e(a) - e(c) + e(j)
        dem2 = Om(s) - e(j) + e(a) 
!       dem2 = Om(t) - e(j) + e(c) 
        dem3 = w - e(a) - Om(t) 
!       dem3 = w - e(c) - Om(s) 
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig + 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z - 2d0*num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              - 2d0*num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)

        dem1 = w - e(a) - e(c) + e(j)
        dem2 = w - e(a) - Om(t)
        dem3 = w - e(c) - Om(s) 
 
        reg1 = (1d0 - exp(-2d0*flow*dem1*dem1))
        reg2 = (1d0 - exp(-2d0*flow*dem2*dem2))
        reg3 = (1d0 - exp(-2d0*flow*dem3*dem3))

        Sig = Sig - num*reg1*reg2*reg3/(dem1*dem2*dem3)

        Z = Z + num*reg1*reg2*reg3/(dem1*dem1*dem2*dem3) &
              + num*reg1*reg2*reg3/(dem1*dem2*dem2*dem3) &
              + num*reg1*reg2*reg3/(dem1*dem2*dem3*dem3)
 
      end do
      end do
    end do
    end do
    end do

  end if
 
end subroutine 
