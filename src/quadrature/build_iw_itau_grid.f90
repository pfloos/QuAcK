
subroutine build_iw_itau_grid(nBas,nOrb,nO,ntimes,nfreqs,verbose,wweight,wcoord,tweight,tcoord,cHF,eHF)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)             :: nBas
  integer,intent(in)             :: verbose
  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nO
  double precision,intent(in)    :: cHF(nBas,nOrb)

! Local variables

  integer                        :: itau
  integer                        :: ifreq
  integer                        :: ibas,jbas
  integer                        :: kind_int

  double precision               :: chem_pot,teval,weval,norm,max_weval
  double precision               :: max_teval_plus,max_teval_minus,max_teval
  double precision               :: alpha,beta,lim_inf,lim_sup

  complex*16                     :: weval_cpx
  complex*16,allocatable         :: G_test(:,:)

! Output variables

  integer,intent(inout)          :: ntimes
  integer,intent(inout)          :: nfreqs
  double precision,intent(inout) :: wweight(nfreqs)
  double precision,intent(inout) :: wcoord(nfreqs)
  double precision,intent(inout) :: tweight(ntimes)
  double precision,intent(inout) :: tcoord(ntimes)
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
   write(*,'(*(a,f20.7))') ' G_test(i tau) ',teval
   do ibas=1,nBas
    write(*,'(*(f15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,f20.7)') ' Largest  tau value for a significant G(+itau) ',teval
 write(*,'(a,f20.7)') '                                 Norm G(+itau) ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,f20.7))') ' G_test(i tau) ',teval
  do ibas=1,nBas
   write(*,'(*(f15.8))') G_test(ibas,:)
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
   write(*,'(*(a,f20.7))') ' G_test(i tau) ',teval
   do ibas=1,nBas
    write(*,'(*(f15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,f20.7)') ' Largest -tau value for a significant G(-itau) ',teval
 write(*,'(a,f20.7)') '                                 Norm G(-itau) ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,f20.7))') ' G_test(i tau) ',teval
  do ibas=1,nBas
   write(*,'(*(f15.8))') G_test(ibas,:)
  enddo
 endif
 max_teval_minus=teval
 
 max_teval=max(abs(max_teval_minus),abs(max_teval_plus))
 write(*,'(a,f20.7)') ' Largest |tau| value for a significant G ',max_teval

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
   write(*,'(*(a,f20.7))') ' G_test(i w) ',weval
   do ibas=1,nBas
    write(*,'(*(f15.8))') G_test(ibas,:)
   enddo
  endif
 enddo
 write(*,'(a,f20.7)') ' Largest w value for a significant G(iw) ',weval
 write(*,'(a,f20.7)') '                           Norm G(-itau) ',norm
 if(verbose/=0) then
  write(*,*)
  write(*,'(*(a,f20.7))') ' G_test(i w) ',weval
  do ibas=1,nBas
   write(*,'(*(f15.8))') G_test(ibas,:)
  enddo
 endif
 max_weval=weval
 write(*,'(a,f20.7)') ' Largest |w| value for a significant G ',max_weval

 
 eHF(:) = eHF(:)+chem_pot
 deallocate(G_test)

!-------------------------!
! Prepare time Quadrature !
!-------------------------!
 kind_int = 1
 lim_inf = 0d0; lim_sup = max_teval;
 alpha = 0d0;  beta  = 0d0;
 call cgqf(ntimes,kind_int,alpha,beta,lim_inf,lim_sup,tcoord,tweight)
! tweight(:)=tweight(:)/((1d0-tcoord(:))**2d0)
! tcoord(:)=tcoord(:)/(1d0-tcoord(:))

!-------------------------------------!
! Prepare second frequency Quadrature !
!-------------------------------------!
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 call cgqf(nfreqs,kind_int,alpha,beta,lim_inf,lim_sup,wcoord,wweight)
 wweight(:)=wweight(:)/((1d0-wcoord(:))**2d0)
 wcoord(:)=wcoord(:)/(1d0-wcoord(:))

end subroutine

