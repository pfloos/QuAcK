
! ---

subroutine  reshufle_hfb(nBas2,nOrb2,c_hfb,eHFB,project) 

! Print one-electron energies and other stuff for G0W0

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas2, nOrb2
  double precision,intent(inout)     :: eHFB(nOrb2)
  double precision,intent(inout)     :: project(nOrb2)
  double precision,intent(inout)     :: c_hfb(nBas2,nOrb2)

! Local variables

  integer                            :: iorb,iorb2,ipos
  double precision                   :: val
  double precision,allocatable       :: eHFB_tmp(:)
  double precision,allocatable       :: project_tmp(:)
  double precision,allocatable       :: c_tmp(:,:)

  allocate(c_tmp(nBas2,nOrb2),project_tmp(nOrb2),eHFB_tmp(nOrb2))

!write(*,*) 'e_HFB'
!do iorb=1,nOrb2
!write(*,*) eHFB(iorb),project(iorb)
!enddo

  ! Reshufle using MOM
  do iorb=1,nOrb2
   val=-1d0
   ipos=0
   do iorb2=1,nOrb2
    if(project(iorb2)>val) then
      val=project(iorb2)
      ipos=iorb2
    endif      
   enddo
   project(ipos)=-1.0d1
   project_tmp(iorb)=val
   eHFB_tmp(iorb)=eHFB(ipos)
   c_tmp(:,iorb)=c_hfb(:,ipos)
  enddo

  ! Reshufle using energies
  do iorb=1,nOrb2/2 ! electronic 1-RDM
   val=1d10
   ipos=0
   do iorb2=1,nOrb2/2
    if(eHFB_tmp(iorb2)<val) then
      val=eHFB_tmp(iorb2)
      ipos=iorb2
    endif      
   enddo
   eHFB_tmp(ipos)=2d10
   eHFB(iorb)=val
   project(iorb)=project_tmp(ipos)
   c_hfb(:,iorb)=c_tmp(:,ipos)
  enddo
  do iorb=nOrb2/2+1,nOrb2 ! hole 1-RDM
   val=1d10
   ipos=0
   do iorb2=nOrb2/2+1,nOrb2
    if(eHFB_tmp(iorb2)<val) then
      val=eHFB_tmp(iorb2)
      ipos=iorb2
    endif      
   enddo
   eHFB_tmp(ipos)=2d10
   eHFB(iorb)=val
   project(iorb)=project_tmp(ipos)
   c_hfb(:,iorb)=c_tmp(:,ipos)
  enddo
  
!write(*,*) 'e_HFB'
!do iorb=1,nOrb2/2
!write(*,*) eHFB(iorb),project(iorb)
!enddo
!write(*,*) '------'
!do iorb=nOrb2/2+1,nOrb2
!write(*,*) eHFB(iorb),project(iorb)
!enddo

 deallocate(c_tmp,project_tmp,eHFB_tmp)

end subroutine 
