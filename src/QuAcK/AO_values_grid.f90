subroutine AO_values_grid(nBas,nShell,CenterShell,TotAngMomShell,KShell,DShell,ExpShell, & 
                          nGrid,root,AO,dAO)

! Compute values of the AOs and their derivatives with respect to the cartesian coordinates

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nShell
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell)
  integer,intent(in)            :: KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK)
  double precision,intent(in)   :: ExpShell(maxShell,maxK)
  double precision,intent(in)   :: root(3,nGrid)
  integer,intent(in)            :: nGrid

! Local variables

  integer                       :: atot,nShellFunction,a(3)
  integer,allocatable           :: ShellFunction(:,:)
  double precision              :: rASq,xA,yA,zA,norm_coeff,prim

  integer                       :: iSh,iShF,iK,iG,iBas

! Output variables

  double precision,intent(out)  :: AO(nBas,nGrid)
  double precision,intent(out)  :: dAO(3,nBas,nGrid)

! Initialization

  iBas = 0
  AO(:,:)    = 0d0
  dAO(:,:,:) = 0d0

!------------------------------------------------------------------------
! Loops over shells
!------------------------------------------------------------------------
  do iSh=1,nShell

    atot = TotAngMomShell(iSh)
    nShellFunction = (atot*atot + 3*atot + 2)/2
    allocate(ShellFunction(1:nShellFunction,1:3))
    call generate_shell(atot,nShellFunction,ShellFunction)

    do iShF=1,nShellFunction

      iBas = iBas + 1
      a(:) = ShellFunction(iShF,:)

      do iG=1,nGrid

        xA = root(1,iG) - CenterShell(iSh,1)
        yA = root(2,iG) - CenterShell(iSh,2)
        zA = root(3,iG) - CenterShell(iSh,3)

!       Calculate distance for exponential

        rASq = xA**2 + yA**2 + zA**2

!------------------------------------------------------------------------
!         Loops over contraction degrees
!-------------------------------------------------------------------------
        do iK=1,KShell(iSh)
          
!         Calculate the exponential part

          prim = DShell(iSh,iK)*norm_coeff(ExpShell(iSh,iK),a)*exp(-ExpShell(iSh,iK)*rASq)        
          AO(iBas,iG) = AO(iBas,iG) + prim

          prim = -2d0*ExpShell(iSh,iK)*prim
          dAO(:,iBas,iG) = dAO(:,iBas,iG) + prim
        
        enddo

        dAO(1,iBas,iG) = xA**(a(1)+1)*yA**a(2)*zA**a(3)*dAO(1,iBas,iG)
        if(a(1) > 0) dAO(1,iBas,iG) = dAO(1,iBas,iG) + dble(a(1))*xA**(a(1)-1)*yA**a(2)*zA**a(3)*AO(iBas,iG)
        
        dAO(2,iBas,iG) = xA**a(1)*yA**(a(2)+1)*zA**a(3)*dAO(2,iBas,iG)
        if(a(2) > 0) dAO(2,iBas,iG) = dAO(2,iBas,iG) + dble(a(2))*xA**a(1)*yA**(a(2)-1)*zA**a(3)*AO(iBas,iG)
        
        dAO(3,iBas,iG) = xA**a(1)*yA**a(2)*zA**(a(3)+1)*dAO(3,iBas,iG)
        if(a(3) > 0) dAO(3,iBas,iG) = dAO(3,iBas,iG) + dble(a(3))*xA**a(1)*yA**a(2)*zA**(a(3)-1)*AO(iBas,iG)

!       Calculate polynmial part

        AO(iBas,iG) = xA**a(1)*yA**a(2)*zA**a(3)*AO(iBas,iG)
        
      enddo
      
    enddo
    deallocate(ShellFunction)
  enddo
!------------------------------------------------------------------------
! End loops over shells 
!------------------------------------------------------------------------

end subroutine AO_values_grid
