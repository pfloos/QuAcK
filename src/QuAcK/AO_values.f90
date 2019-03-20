subroutine AO_values(doDrift,nBas,nShell,nWalk,CenterShell,TotAngMomShell,KShell,DShell,ExpShell,r,AO,dAO)

! Compute values of the AOs and their derivatives (if required)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doDrift
  integer,intent(in)            :: nBas,nShell,nWalk
  double precision,intent(in)   :: CenterShell(maxShell,3)
  integer,intent(in)            :: TotAngMomShell(maxShell),KShell(maxShell)
  double precision,intent(in)   :: DShell(maxShell,maxK),ExpShell(maxShell,maxK)
  double precision,intent(in)   :: r(nWalk,2,3)

! Local variables

  integer                       :: atot,nShellFunction,a(3)
  integer,allocatable           :: ShellFunction(:,:)
  double precision              :: rASq,xA,yA,zA,norm_coeff,prim

  integer                       :: iSh,iShF,iK,iW,iEl,iBas,ixyz

! Output variables

  double precision,intent(out)  :: AO(nWalk,2,nBas),dAO(nWalk,2,3,nBas)

! Initialization

  AO = 0d0
  if(doDrift) dAO = 0d0
  iBas = 0

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
      a(1) = ShellFunction(iShF,1)
      a(2) = ShellFunction(iShF,2)
      a(3) = ShellFunction(iShF,3)

      do iW=1,nWalk
        do iEl=1,2

          xA = r(iW,iEl,1) - CenterShell(iSh,1)
          yA = r(iW,iEl,2) - CenterShell(iSh,2)
          zA = r(iW,iEl,3) - CenterShell(iSh,3)

!         Calculate distance for exponential

          rASq = xA**2 + yA**2 + zA**2

!------------------------------------------------------------------------
!         Loops over contraction degrees
!-------------------------------------------------------------------------
          do iK=1,KShell(iSh)
            
!           Calculate the exponential part
            prim = DShell(iSh,iK)*norm_coeff(ExpShell(iSh,iK),a)*exp(-ExpShell(iSh,iK)*rASq)        
            AO(iW,iEl,iBas) = AO(iW,iEl,iBas) + prim

            if(doDrift) then 
              prim = -2d0*ExpShell(iSh,iK)*prim
              do ixyz=1,3
                dAO(iW,iEl,ixyz,iBas) = dAO(iW,iEl,ixyz,iBas) + prim
              enddo
            endif
          
          enddo

          if(doDrift) then         

            dAO(iW,iEl,1,iBas) = xA**(a(1)+1)*yA**a(2)*zA**a(3)*dAO(iW,iEl,1,iBas)
            if(a(1) > 0) dAO(iW,iEl,1,iBas) = dAO(iW,iEl,1,iBas) + dble(a(1))*xA**(a(1)-1)*yA**a(2)*zA**a(3)*AO(iW,iEl,iBas)
            
            dAO(iW,iEl,2,iBas) = xA**a(1)*yA**(a(2)+1)*zA**a(3)*dAO(iW,iEl,2,iBas)
            if(a(2) > 0) dAO(iW,iEl,2,iBas) = dAO(iW,iEl,2,iBas) + dble(a(2))*xA**a(1)*yA**(a(2)-1)*zA**a(3)*AO(iW,iEl,iBas)
            
            dAO(iW,iEl,3,iBas) = xA**a(1)*yA**a(2)*zA**(a(3)+1)*dAO(iW,iEl,3,iBas)
            if(a(3) > 0) dAO(iW,iEl,3,iBas) = dAO(iW,iEl,3,iBas) + dble(a(3))*xA**a(1)*yA**a(2)*zA**(a(3)-1)*AO(iW,iEl,iBas)
            
          endif

!         Calculate polynmial part

          AO(iW,iEl,iBas) = xA**a(1)*yA**a(2)*zA**a(3)*AO(iW,iEl,iBas)

        enddo
      enddo
      
    enddo
    deallocate(ShellFunction)
  enddo
!------------------------------------------------------------------------
! End loops over shells 
!------------------------------------------------------------------------

end subroutine AO_values
