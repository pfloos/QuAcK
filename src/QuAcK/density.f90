subroutine density(doDrift,nBas,nWalk,P,gAO,dgAO,g,dg)

! Calculate the Green functions

  implicit none

! Input variables

  logical,intent(in)            :: doDrift
  integer,intent(in)            :: nBas,nWalk
  double precision,intent(in)   :: P(nBas,nBas),gAO(nWalk,2,nBas),dgAO(nWalk,2,3,nBas)

! Local variables

  integer                       :: iW,iEl,ixyz,mu,nu

! Output variables

  double precision,intent(out)  :: g(nWalk,2),dg(nWalk,2,3)

  g = 0d0
  do iW=1,nWalk
   do iEl=1,2
     do mu=1,nBas
       do nu=1,nBas
         g(iW,iEl) = g(iW,iEl) + gAO(iW,iEl,mu)*P(mu,nu)*gAO(iW,iEl,nu)
       enddo
     enddo
   enddo
  enddo

  if(doDrift) then

    dg = 0d0
    do iW=1,nWalk
      do iEl=1,2
        do ixyz=1,3
          do mu=1,nBas
            do nu=1,nBas
              dg(iW,iEl,ixyz) = dg(iW,iEl,ixyz)                                & 
                              + P(mu,nu)*(dgAO(iW,iEl,ixyz,mu)*gAO(iW,iEl,nu)  &
                              +            gAO(iW,iEl,mu)*dgAO(iW,iEl,ixyz,nu))
            enddo
          enddo
        enddo
      enddo
    enddo

  endif

end subroutine density
