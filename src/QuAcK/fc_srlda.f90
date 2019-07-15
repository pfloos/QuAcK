subroutine fc_srlda(nBas,nGrid,weight,MO,rho,mu,eG0W0)

! Compute sr-lda ec

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: MO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: mu(nGrid)
  double precision,intent(in)   :: eG0W0(nBas)

! Local variables

  integer                       :: iG,p
  double precision              :: r,ra,rb,rap,ram
  double precision              :: rs,rsp,rsm
  double precision              :: ec,ecp,ecm,vcup,vcdw
  double precision,parameter    :: delta = 1d-6
  double precision,allocatable  :: de(:)

! Memory allocation

  allocate(de(nBas))

! Initialization

  de(:) = 0d0

  do iG=1,ngrid

    r  = max(0d0,rho(iG))
    ra = 0.5d0*r
    rb = 0.5d0*r

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

!     call lsdsr(rs,0d0,mu(iG),ec,vcup,vcdw)
      if(abs(ra) > delta) then 

        rap = ra + delta
        ram = ra - delta

        rsp = (4d0*pi*rap/3d0)**(-1d0/3d0)
        rsm = (4d0*pi*ram/3d0)**(-1d0/3d0)

!      call lsdsr(rsp,0d0,mu(iG),ecp,vcup,vcdw)
!      call lsdsr(rsm,0d0,mu(iG),ecm,vcup,vcdw)
       call lsdsr(rs,0d0,mu(iG),ec,vcup,vcdw)


!      call ESRC_MD_LDAERF(mu(iG),rap,rb,.true.,ecp)
!      call ESRC_MD_LDAERF(mu(iG),ram,rb,.true.,ecm)

!       vcup = (ecp - ecm)/(2d0*delta)

      else

        vcup = 0d0  

      end if

      do p=1,nBas

        de(p)= de(p) + weight(iG)*vcup*MO(p,iG)**2

      end do
 
    end if

  end do

  print*, 'Eigenvalues correction from srDFT (in eV)'
  call matout(nBas,1,de(:)*HaToeV)

  print*, 'Corrected G0W0 eigenvalues (in eV)'
  call matout(nBas,1,(eG0W0(:) + de(:))*HaToeV)

  

end subroutine fc_srlda
