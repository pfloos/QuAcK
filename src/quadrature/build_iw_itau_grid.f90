
subroutine build_iw_itau_grid(nBas,nOrb,nO,ntimes,nfreqs,verbose,cHF,eHF)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)             :: nBas
  integer,intent(in)             :: verbose
  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nO
  double precision,intent(in)    :: cHF(nBas,nOrb)

! Local variables

  integer                        :: iunit=311
  integer                        :: itau
  integer                        :: ifreq
  integer                        :: ibas,jbas
  integer                        :: kind_int
  integer                        :: ntimes_01,ntimes_intervals

  double precision               :: m,y0
  double precision               :: chem_pot,teval,weval,norm,max_weval
  double precision               :: max_teval_plus,max_teval_minus,max_teval
  double precision               :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable   :: tcoord_01(:),tweight_01(:)
  double precision,allocatable   :: wweight(:),wcoord(:)
  double precision,allocatable   :: tweight(:),tcoord(:)

  complex*16                     :: weval_cpx
  complex*16,allocatable         :: G_test(:,:)

! Output variables

  integer,intent(inout)          :: ntimes
  integer,intent(inout)          :: nfreqs
  double precision,intent(inout) :: eHF(nOrb)

!------------------------------------------------------------------------
! Build iw and itau grids 
!------------------------------------------------------------------------

 allocate(G_test(nBas,nBas))
 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 eHF(:) = eHF(:)-chem_pot

 ! Find the largest +tau for Go(i tau)
 teval=2d0
 do
  call G0itau_ao_RHF(nBas,nOrb,nO,teval,G_test,cHF,eHF)
  norm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    norm=norm+abs(G_test(ibas,jbas))
   enddo
  enddo
  if(norm<1d-8) exit
  teval=teval+2d0
  if(verbose/=0) then
   write(*,*)
   write(*,'(*(a,e20.7))') ' G_test(i tau) ',teval
   do ibas=1,nBas
    write(*,'(*(e15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,e20.7)') ' Largest  tau value for a significant G(+itau) ',teval
 write(*,'(a,e20.7)') '                                 Norm G(+itau) ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,e20.7))') ' G_test(i tau) ',teval
  do ibas=1,nBas
   write(*,'(*(e15.8))') G_test(ibas,:)
  enddo
 endif
 max_teval_plus=teval

 ! Find the largest -tau for Go(-i tau)
 teval=-2d0
 do
  call G0itau_ao_RHF(nBas,nOrb,nO,teval,G_test,cHF,eHF)
  norm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    norm=norm+abs(G_test(ibas,jbas))
   enddo
  enddo
  if(norm<1d-8) exit
  teval=teval-2d0
  if(verbose/=0) then
   write(*,*)
   write(*,'(*(a,e20.7))') ' G_test(i tau) ',teval
   do ibas=1,nBas
    write(*,'(*(e15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,e20.7)') ' Largest -tau value for a significant G(-itau) ',teval
 write(*,'(a,e20.7)') '                                 Norm G(-itau) ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,e20.7))') ' G_test(i tau) ',teval
  do ibas=1,nBas
   write(*,'(*(e15.8))') G_test(ibas,:)
  enddo
 endif
 max_teval_minus=teval
 
 max_teval=max(abs(max_teval_minus),abs(max_teval_plus))
 write(*,'(a,e20.7)') ' Largest |tau| value for a significant G ',max_teval

 ! Find the largest w for Go(i w)
 weval=1d4
 do
  weval_cpx=weval*im
  call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_test)
  norm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    norm=norm+abs(G_test(ibas,jbas))
   enddo
  enddo
  if(norm<1d-6) exit
  weval=weval+1d4
  if(verbose/=0) then
   write(*,*)
   write(*,'(*(a,e20.7))') ' G_test(i w) ',weval
   do ibas=1,nBas
    write(*,'(*(e15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,e20.7)') ' Largest   w  value for a significant  G(+iw)  ',weval
 write(*,'(a,e20.7)') '                                  Norm G(+iw)  ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,e20.7))') ' G_test(i w) ',weval
  do ibas=1,nBas
   write(*,'(*(e15.8))') G_test(ibas,:)
  enddo
 endif
 max_weval=weval
 write(*,'(a,e20.7)') ' Largest  |w|  value for a significant G ',max_weval

 
 eHF(:) = eHF(:)+chem_pot
 deallocate(G_test)

!-------------------------!
! Prepare time Quadrature !
!-------------------------!
! ntimes_01 = 50 ! From 0 to 1 we take 50 points always
 ntimes_01 =1000 ! TEST TODO
 ntimes_intervals = 1
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1;
 alpha = 0d0;  beta  = 0d0;
 allocate(tweight_01(ntimes_01),tcoord_01(ntimes_01))
 call cgqf(ntimes_01,kind_int,alpha,beta,lim_inf,lim_sup,tcoord_01,tweight_01)

 ntimes=ntimes_01*ntimes_intervals
 allocate(tweight(ntimes),tcoord(ntimes))
 
 m=max_teval;y0=0d0;
 tweight(:)=m*tweight_01(:)
 tcoord(:)=m*tcoord_01(:)+y0

 open(unit=iunit, form='formatted', file='tcoord.txt')
 do itau=1,ntimes
   write(iunit,'(f50.15)') tcoord(itau)
 enddo 
 close(iunit)
 open(unit=iunit, form='formatted', file='tweight.txt')
 do itau=1,ntimes
   write(iunit,'(f50.15)') tweight(itau)
 enddo 
 close(iunit)
 deallocate(tweight_01,tcoord_01)
 deallocate(tweight,tcoord)

!------------------------------!
! Prepare frequency Quadrature !
!------------------------------!

 ! TODO
 nfreqs=800
 allocate(wweight(nfreqs),wcoord(nfreqs))


 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 call cgqf(nfreqs,kind_int,alpha,beta,lim_inf,lim_sup,wcoord,wweight)
 wweight(:)=wweight(:)/((1d0-wcoord(:))**2d0)
 wcoord(:)=wcoord(:)/(1d0-wcoord(:))

 open(unit=iunit, form='formatted', file='wcoord.txt')
 do ifreq=1,nfreqs
   write(iunit,'(f50.15)') wcoord(ifreq)
 enddo 
 close(iunit)
 open(unit=iunit, form='formatted', file='wweight.txt')
 do ifreq=1,nfreqs
   write(iunit,'(f50.15)') wweight(ifreq)
 enddo 
 close(iunit)
 deallocate(wweight,wcoord)

end subroutine

