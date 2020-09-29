subroutine transition_probability(nWalk,dt,D,r,rp,F,Fp,T,Tp)

! Compute transition probability

  implicit none

! Input variables

  integer,intent(in)            :: nWalk
  double precision,intent(in)   :: dt,D
  double precision,intent(in)   :: r(nWalk,1:2,1:3), F(nWalk,1:2,1:3)
  double precision,intent(in)   :: rp(nWalk,1:2,1:3),Fp(nWalk,1:2,1:3)

! Local variables 
  
  integer               :: iW,iEl,ixyz

! Output variables

  double precision,intent(out)  :: T(nWalk),Tp(nWalk)

! Initialize 

  T  = 0d0
  Tp = 0d0

! Compute 

  do iW=1,nWalk
    do iEl=1,2
      do ixyz=1,3
        T(iW)  = T(iW)  + (rp(iW,iEl,ixyz) - r(iW,iEl,ixyz)  - D*dt*F(iW,iEl,ixyz))**2
        Tp(iW) = Tp(iW) + (r(iW,iEl,ixyz)  - rp(iW,iEl,ixyz) - D*dt*Fp(iW,iEl,ixyz))**2
      enddo
    enddo
  enddo

  T(:)  = exp(-0.25d0*T(:)/(D*dt))
  Tp(:) = exp(-0.25d0*Tp(:)/(D*dt))

end subroutine transition_probability
