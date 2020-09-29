subroutine NDrift(nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,P,r,w,F)

! Compute quantum force numerically

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nWalk,nBas
  double precision,intent(in)   :: P(nBas,nBas),r(nWalk,2,3),w(nWalk)

  integer,intent(in)            :: nShell
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: CenterShell(maxShell,3),DShell(maxShell,maxK),ExpShell(maxShell,maxK)

! Local variables

  integer                       :: iW,iEl,ixyz
  double precision              :: delta
  double precision              :: wp,wm
  double precision,allocatable  :: rp(:,:,:),rm(:,:,:),r12p(:),r12m(:)
  double precision,allocatable  :: gAOp(:,:,:),dgAOp(:,:,:,:),gAOm(:,:,:),dgAOm(:,:,:,:)
  double precision,allocatable  :: gp(:,:),dgp(:,:,:),gm(:,:),dgm(:,:,:)


! Output variables

  double precision,intent(out)  :: F(nWalk,2,3)

  allocate(rp(nWalk,2,3),rm(nWalk,2,3),r12p(nWalk),r12m(nWalk), &
           gAOp(nWalk,2,nBas),dgAOp(nWalk,2,3,nBas),gAOm(nWalk,2,nBas),dgAOm(nWalk,2,3,nBas), &
           gp(nWalk,2),dgp(nWalk,2,3),gm(nWalk,2),dgm(nWalk,2,3))

  do iW=1,nWalk
    do iEl=1,2
      do ixyz=1,3

        delta = 1d-6

        rp = r
        rm = r

        rp(iW,iEl,ixyz) = r(iW,iEl,ixyz) + delta
        rm(iW,iEl,ixyz) = r(iW,iEl,ixyz) - delta
        
        call AO_values(.false.,nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,rp,gAOp,dgAOp)
        call AO_values(.false.,nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,rm,gAOm,dgAOm)
        
        call Density(.false.,nBas,nWalk,P,gAOp,dgAOp,gp,dgp)
        call Density(.false.,nBas,nWalk,P,gAOm,dgAOm,gm,dgm)
        
        call rij(nWalk,rp,r12p)
        call rij(nWalk,rm,r12m)
        
        wp = gp(iW,1)*gp(iW,2)/r12p(iW)
        wm = gm(iW,1)*gm(iW,2)/r12m(iW)
        
        F(iW,iEl,ixyz) = (wp - wm)/(2d0*delta*w(iw))
      enddo
    enddo
  enddo

! print*,'NF',F


end subroutine NDrift
