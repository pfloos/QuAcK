
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

  logical                        :: all_deactivated
  logical,allocatable            :: interval_todo(:)
  logical,allocatable            :: interval_todo2(:)

  integer                        :: iunit=311
  integer                        :: iter,igrid,jgrid,iinterval,jinterval
  integer                        :: itau
  integer                        :: ifreq
  integer                        :: ibas,jbas
  integer                        :: kind_int
  integer                        :: n01,ntimes_t,nintervals,nintervals_old,nintervals_twice,ngrid,niter_max

  double precision               :: ran_num
  double precision               :: thrs_tnorm,thrs_wnorm
  double precision               :: thrs_interval,max_err_G_set,err_G_set,err_G_set_old
  double precision               :: chem_pot,teval,weval,tweight_eval,wweight_eval,norm,max_weval
  double precision               :: max_teval_plus,max_teval_minus,max_teval,min_teval
  double precision               :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable   :: coord_01(:),weight_01(:)
  double precision,allocatable   :: wweight(:),wcoord(:)
  double precision,allocatable   :: tweight(:),tcoord(:)
  double precision,allocatable   :: interval_r(:,:)
  double precision,allocatable   :: interval_tmp_weight(:,:),interval_tmp_coord(:,:)
  double precision,allocatable   :: interval_r2(:,:)
  double precision,allocatable   :: interval_tmp_weight2(:,:),interval_tmp_coord2(:,:)

  complex*16                     :: weval_cpx
  complex*16,allocatable         :: G_test(:,:),G_tmp1(:,:)
  complex*16,allocatable         :: G_set_test(:,:,:)
  complex*16,allocatable         :: interval_vals(:,:,:)

! Output variables

  integer,intent(inout)          :: ntimes
  integer,intent(inout)          :: nfreqs
  double precision,intent(inout) :: eHF(nOrb)

!------------------------------------------------------------------------
! Build iw and itau grids 
!------------------------------------------------------------------------
 
 thrs_tnorm=1d-6
 thrs_wnorm=1d-6
 min_teval=2d-4
 ngrid=40
 allocate(G_set_test(ngrid,nBas,nBas))
 allocate(G_test(nBas,nBas),G_tmp1(nbas,nBas))
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
  if(norm<thrs_tnorm) exit
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
  if(norm<thrs_tnorm) exit
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
  if(norm<thrs_wnorm) exit
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
 write(*,*) 

!-------------------------!
! Prepare time Quadrature !
!-------------------------!
 niter_max = 40
 n01 = 50 ! From 0 to 1 we always take 50 points
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 allocate(weight_01(n01),coord_01(n01))
 call cgqf(n01,kind_int,alpha,beta,lim_inf,lim_sup,coord_01,weight_01)
 ! Build 40 trial freqs from 0 to max_weval
 allocate(wcoord(ngrid)) 
 wcoord(1)=0d0
 wcoord(ngrid)=max_weval
 do ifreq=2,ngrid-1
  call random_number(ran_num)
  ran_num=ran_num**1d1
  ran_num=ran_num*max_weval
  wcoord(ifreq)=ran_num 
 enddo
 call sort_ascending(ngrid,wcoord)
 ! Prepare the reference G(iw) values
 do ifreq=1,ngrid
  weval_cpx=wcoord(ifreq)*im
  call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_test)
  G_set_test(ifreq,:,:)=G_test(:,:)
 enddo
 ! Optimize the weights and coordinates with an adaptative quadrature
 iter=-1
 nintervals=1
 max_err_G_set=5d-4
 thrs_interval=1d-3
 allocate(tweight(1),tcoord(1))
 allocate(interval_vals(ngrid,nBas,nBas))
 allocate(interval_tmp_weight(nintervals,n01))
 allocate(interval_tmp_coord(nintervals,n01))
 allocate(interval_r(nintervals,2))
 allocate(interval_todo(nintervals))
 iinterval=1
 interval_todo(:)=.true.
 interval_r(iinterval,1)=min_teval ! Start tau grid at least at 1d-4 because when tau->0 is not continuos [ Lim tau->0+ /= Lim tau->0- ]
 interval_r(iinterval,2)=max_teval
 interval_tmp_weight(iinterval,:)=weight_01(:)*(interval_r(iinterval,2)-interval_r(iinterval,1))
 interval_tmp_coord(iinterval,:) = coord_01(:)*(interval_r(iinterval,2)-interval_r(iinterval,1))+interval_r(iinterval,1)
 do
  iter=iter+1

  ! Update intervals
  if(iter>0) then
   nintervals_twice=2*nintervals
   allocate(interval_tmp_weight2(nintervals_twice,n01))
   allocate(interval_tmp_coord2(nintervals_twice,n01))
   allocate(interval_r2(nintervals_twice,2))
   allocate(interval_todo2(nintervals_twice))
   interval_todo2(:)=.true.
   do iinterval=1,nintervals
    if(interval_todo(iinterval)) then
     ! For the ith active interval [a,b] compute all G_int(iw_k) = Conjg[ -i int _a ^b [ G(i tau) Exp(i tau w_k)  + G(-i tau) Exp(-i tau w_k) ] d tau ]
     interval_vals=czero 
     do ifreq=1,ngrid
      G_test=czero
      do itau=1,n01
       teval=interval_tmp_coord(iinterval,itau)
       tweight_eval=interval_tmp_weight(iinterval,itau)
       call G0itau_ao_RHF(nBas,nOrb,nO, teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp( im*teval*wcoord(ifreq))  ! G(i tau)
       call G0itau_ao_RHF(nBas,nOrb,nO,-teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp(-im*teval*wcoord(ifreq))  ! G(-i tau)
      enddo
      G_test=conjg(G_test)
      interval_vals(ifreq,:,:)=G_test(:,:)
     enddo
     ! Divide the ith interval in halves and compute coords and weights
     interval_r2(2*iinterval-1,1)=interval_r(iinterval,1) 
     interval_r2(2*iinterval-1,2)=interval_r(iinterval,1)+0.5d0*(interval_r(iinterval,2)-interval_r(iinterval,1)) 
     interval_r2(2*iinterval  ,1)=interval_r2(2*iinterval-1,2) 
     interval_r2(2*iinterval  ,2)=interval_r(iinterval,2)
     interval_tmp_weight2(2*iinterval-1,:)=weight_01(:)*(interval_r2(2*iinterval-1,2)-interval_r2(2*iinterval-1,1))
     interval_tmp_coord2(2*iinterval-1,:) = coord_01(:)*(interval_r2(2*iinterval-1,2)-interval_r2(2*iinterval-1,1))+interval_r2(2*iinterval-1,1)
     interval_tmp_weight2(2*iinterval  ,:)=weight_01(:)*(interval_r2(2*iinterval  ,2)-interval_r2(2*iinterval  ,1))
     interval_tmp_coord2(2*iinterval  ,:) = coord_01(:)*(interval_r2(2*iinterval  ,2)-interval_r2(2*iinterval  ,1))+interval_r2(2*iinterval  ,1)
     ! Check the difference in G_int(iw_k) computed as [a,b] vs [a,0.5(b-a)+a] + [0.5(b-a)+a,b]
     err_G_set=0d0
     do ifreq=1,ngrid
      G_test=czero
      do itau=1,n01
       ! [a,0.5(b-a)+a]
       teval=interval_tmp_coord2(2*iinterval-1,itau)
       tweight_eval=interval_tmp_weight2(2*iinterval-1,itau)
       call G0itau_ao_RHF(nBas,nOrb,nO, teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp( im*teval*wcoord(ifreq))  ! G(i tau)
       call G0itau_ao_RHF(nBas,nOrb,nO,-teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp(-im*teval*wcoord(ifreq))  ! G(-i tau)
       ! [0.5(b-a)+a,b]
       teval=interval_tmp_coord2(2*iinterval  ,itau)
       tweight_eval=interval_tmp_weight2(2*iinterval  ,itau)
       call G0itau_ao_RHF(nBas,nOrb,nO, teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp( im*teval*wcoord(ifreq))  ! G(i tau)
       call G0itau_ao_RHF(nBas,nOrb,nO,-teval,G_tmp1,cHF,eHF)
       G_test(:,:)=G_test(:,:)-im*tweight_eval*G_tmp1(:,:)*Exp(-im*teval*wcoord(ifreq))  ! G(-i tau)
      enddo
      G_test=conjg(G_test)
      G_test(:,:)=G_test(:,:)-interval_vals(ifreq,:,:)
      do ibas=1,nBas
       do jbas=1,nBas
        err_G_set=err_G_set+abs(G_test(ibas,jbas))
       enddo
      enddo
     enddo
     ! If we entered here, we always retain the new 2 halves of the original interval
     ! If diff. is small we keep the small intervals, but deactivate these two intervals
     err_G_set=err_G_set/(nBas*nBas*ngrid)
     if(err_G_set<thrs_interval) then
      interval_todo2(2*iinterval-1)=.false.
      interval_todo2(2*iinterval  )=.false.
     endif
    else ! The interval is inactive (we store it at the odd row and set to 0 the even row)
     interval_tmp_weight2(2*iinterval-1,:)=interval_tmp_weight(iinterval,:)
     interval_tmp_coord2(2*iinterval-1,:)=interval_tmp_coord(iinterval,:)
     interval_r2(2*iinterval-1,:)=interval_r(iinterval,:)
     interval_todo2(2*iinterval-1)=interval_todo(iinterval)
     interval_tmp_weight2(2*iinterval ,:)=0d0
     interval_tmp_coord2(2*iinterval  ,:)=0d0
     interval_r2(2*iinterval ,:)=0d0
     interval_todo2(2*iinterval)=.false.
    endif 
   enddo
   ! Find the intervals that survived
   nintervals=0
   do iinterval=1,nintervals_twice
    if( abs(interval_r2(iinterval,1))<1d-8 .and. abs(interval_r2(iinterval,2))<1d-8 ) then
     cycle
    else
     nintervals=nintervals+1 
    endif
   enddo
   ! Save in interval_tmp_weight, interval_tmp_coord, interval_r, and interval_todo (i.e., the intervals that survived)
   deallocate(interval_tmp_weight) 
   deallocate(interval_tmp_coord) 
   deallocate(interval_r)
   deallocate(interval_todo)
   allocate(interval_tmp_weight(nintervals,n01))
   allocate(interval_tmp_coord(nintervals,n01))
   allocate(interval_r(nintervals,2))
   allocate(interval_todo(nintervals))
   jinterval=1
   do iinterval=1,nintervals_twice
    if( abs(interval_r2(iinterval,1))<1d-8 .and. abs(interval_r2(iinterval,2))<1d-8 ) then
     cycle
    else
     interval_tmp_weight(jinterval,:)=interval_tmp_weight2(iinterval,:) 
     interval_tmp_coord(jinterval,:)=interval_tmp_coord2(iinterval,:) 
     interval_r(jinterval,:)=interval_r2(iinterval,:) 
     interval_todo(jinterval)=interval_todo2(iinterval)
     jinterval=jinterval+1 
    endif
   enddo
   deallocate(interval_tmp_weight2) 
   deallocate(interval_tmp_coord2) 
   deallocate(interval_r2)
   deallocate(interval_todo2)
  endif

  ! Transfer from the matrices to the vectors tweight and tcoord to compute the full
  !    G_int(iw_k) = Conjg[ -i int _0 ^Infty [ G(i tau) Exp(i tau w_k)  + G(-i tau) Exp(-i tau w_k) ] d tau ]
  ntimes=n01*nintervals
  deallocate(tweight,tcoord)
  allocate(tweight(ntimes),tcoord(ntimes))
  jgrid=1
  do iinterval=1,nintervals
   do igrid=1,n01
    tweight(jgrid) = interval_tmp_weight(iinterval,igrid)
    tcoord(jgrid)  =  interval_tmp_coord(iinterval,igrid)
    jgrid=jgrid+1
   enddo
  enddo

  ! Check the global error for G(i tau) -> G(i w)
  err_G_set=0d0
  do ifreq=1,ngrid
   G_test=czero
   do itau=1,ntimes
    call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_tmp1,cHF,eHF)
    G_test(:,:)=G_test(:,:)-im*tweight(itau)*G_tmp1(:,:)*Exp( im*tcoord(itau)*wcoord(ifreq))  ! G(i tau)
    call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_tmp1,cHF,eHF)
    G_test(:,:)=G_test(:,:)-im*tweight(itau)*G_tmp1(:,:)*Exp(-im*tcoord(itau)*wcoord(ifreq))  ! G(-i tau)
   enddo
   G_test=conjg(G_test)
   if(verbose/=0) then ! Print for comparison and debug
    write(*,'(a,f25.15)') ' wcoord_test',wcoord(ifreq)
    write(*,*) 'Integrated'
    do ibas=1,nbas
     write(*,'(*(f15.8))') G_test(ibas,:)
    enddo
    write(*,*) 'Analytic'
    do ibas=1,nbas
     write(*,'(*(f15.8))') G_set_test(ifreq,ibas,:)
    enddo
    write(*,*)
   endif
   G_test(:,:)=G_test(:,:)-G_set_test(ifreq,:,:)
   do ibas=1,nBas
    do jbas=1,nBas
     err_G_set=err_G_set+abs(G_test(ibas,jbas))
    enddo
   enddo
  enddo
  err_G_set=err_G_set/(nBas*nBas*ngrid) ! MAE
  write(*,'(a,f15.8,a,i5,a,i8,a)') ' MAE in reproducing G(iw) [Integrated - Analytic]',err_G_set,' at iter ',iter,' with ',ntimes,' grid'

  if(iter>0) then
   if(abs(err_G_set_old-err_G_set)<1d-8 .or. nintervals_old==nintervals) then
    write(*,'(a,e15.8,a,e15.8)') ' Reducing the MAE accepted for each interval integrations form ',thrs_interval,' to ',1d-1*thrs_interval 
    thrs_interval=1d-1*thrs_interval
    interval_todo(:)=.true.
   endif
  endif
  err_G_set_old=err_G_set
  nintervals_old=nintervals

  all_deactivated=.true.
  do iinterval=1,nintervals
   if(interval_todo(iinterval)) all_deactivated=.false. 
  enddo
  if(all_deactivated) then
   write(*,'(a,e15.8,a,e15.8)') ' Reducing the MAE accepted for each interval integrations form ',thrs_interval,' to ',1d-1*thrs_interval 
   thrs_interval=1d-1*thrs_interval
   interval_todo(:)=.true.
  endif

  if(err_G_set<max_err_G_set) exit

  if(iter==niter_max) exit

 enddo
 write(*,*) 

 deallocate(interval_vals)
 deallocate(interval_tmp_weight)
 deallocate(interval_tmp_coord)
 deallocate(interval_r,interval_todo)
 deallocate(wcoord)
 
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
 deallocate(weight_01,coord_01)
 deallocate(tweight,tcoord)

!------------------------------!
! Prepare frequency Quadrature !
!------------------------------!
 niter_max = 18
 n01 = 2000 ! From 0 to 1 we always take 2000 points
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 allocate(weight_01(n01),coord_01(n01))
 call cgqf(n01,kind_int,alpha,beta,lim_inf,lim_sup,coord_01,weight_01)
 ! Build 40 trial freqs from -max_teval to max_teval
 allocate(tcoord(ngrid))
 tcoord(1)=-max_teval
 tcoord(2)= max_teval
 tcoord(3)=-min_teval
 tcoord(4)= min_teval
 do itau=2,(ngrid-2)/2
  call random_number(ran_num)
  ran_num=ran_num**5d0
  ran_num=ran_num*(max_teval-min_teval)+min_teval ! Min tau used is 1d-4
  tcoord(2*(itau+1)-1)= ran_num
  tcoord(2*(itau+1)  )=-ran_num
 enddo
 call sort_ascending(ngrid,tcoord)
 ! Prepare the reference G(itau) values
 do itau=1,ngrid
  call G0itau_ao_RHF(nBas,nOrb,nO,tcoord(itau),G_test,cHF,eHF)
  G_set_test(itau,:,:)=G_test(:,:)
 enddo
 ! Optimize the weights and coordinates with an adaptative quadrature
 iter=-1
 nintervals=1
 max_err_G_set=5d-4
 thrs_interval=1d-2
 allocate(wweight(1),wcoord(1))
 allocate(interval_vals(ngrid,nBas,nBas))
 allocate(interval_tmp_weight(nintervals,n01))
 allocate(interval_tmp_coord(nintervals,n01))
 allocate(interval_r(nintervals,2))
 allocate(interval_todo(nintervals))
 iinterval=1
 interval_todo(:)=.true.
 interval_r(iinterval,1)=0d0
 interval_r(iinterval,2)=max_weval
 interval_tmp_weight(iinterval,:)=weight_01(:)*(interval_r(iinterval,2)-interval_r(iinterval,1))
 interval_tmp_coord(iinterval,:) = coord_01(:)*(interval_r(iinterval,2)-interval_r(iinterval,1))+interval_r(iinterval,1)
 do
  iter=iter+1

  ! Update intervals
  if(iter>0) then
   nintervals_twice=2*nintervals
   allocate(interval_tmp_weight2(nintervals_twice,n01))
   allocate(interval_tmp_coord2(nintervals_twice,n01))
   allocate(interval_r2(nintervals_twice,2))
   allocate(interval_todo2(nintervals_twice))
   interval_todo2(:)=.true.
   do iinterval=1,nintervals
    if(interval_todo(iinterval)) then
     ! For the ith active interval [a,b] compute all G_int(itau_k) = i/(2 pi) int _a ^b [ G(i w) Exp(i tau_k w)  + G(-i w) Exp(-i tau_k w) ] d w
     interval_vals=czero 
     do itau=1,ngrid
      G_test=czero
      do ifreq=1,n01
       weval=interval_tmp_coord(iinterval,ifreq)
       wweight_eval=interval_tmp_weight(iinterval,ifreq)
       weval_cpx=weval*im
       call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_tmp1)
       G_tmp1(:,:)=G_tmp1(:,:)*Exp(im*tcoord(itau)*weval)
       G_test(:,:)=G_test(:,:)+im*wweight_eval*(G_tmp1(:,:)+conjg(G_tmp1(:,:))) ! G(i w) + G(-i w)
      enddo
      G_test=G_test/(2d0*pi)
      interval_vals(itau,:,:)=G_test(:,:)
     enddo
     ! Divide the ith interval in halves and compute coords and weights
     interval_r2(2*iinterval-1,1)=interval_r(iinterval,1) 
     interval_r2(2*iinterval-1,2)=interval_r(iinterval,1)+0.5d0*(interval_r(iinterval,2)-interval_r(iinterval,1)) 
     interval_r2(2*iinterval  ,1)=interval_r2(2*iinterval-1,2) 
     interval_r2(2*iinterval  ,2)=interval_r(iinterval,2)
     interval_tmp_weight2(2*iinterval-1,:)=weight_01(:)*(interval_r2(2*iinterval-1,2)-interval_r2(2*iinterval-1,1))
     interval_tmp_coord2(2*iinterval-1,:) = coord_01(:)*(interval_r2(2*iinterval-1,2)-interval_r2(2*iinterval-1,1))+interval_r2(2*iinterval-1,1)
     interval_tmp_weight2(2*iinterval  ,:)=weight_01(:)*(interval_r2(2*iinterval  ,2)-interval_r2(2*iinterval  ,1))
     interval_tmp_coord2(2*iinterval  ,:) = coord_01(:)*(interval_r2(2*iinterval  ,2)-interval_r2(2*iinterval  ,1))+interval_r2(2*iinterval  ,1)
     ! Check the difference in G_int(itau_k) computed as [a,b] vs [a,0.5(b-a)+a] + [0.5(b-a)+a,b]
     err_G_set=0d0
     do itau=1,ngrid
      G_test=czero
      do ifreq=1,n01
       ! [a,0.5(b-a)+a]
       weval=interval_tmp_coord2(2*iinterval-1,ifreq)
       wweight_eval=interval_tmp_weight2(2*iinterval-1,ifreq)
       weval_cpx=weval*im
       call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_tmp1)
       G_tmp1(:,:)=G_tmp1(:,:)*Exp(im*tcoord(itau)*weval)
       G_test(:,:)=G_test(:,:)+im*wweight_eval*(G_tmp1(:,:)+conjg(G_tmp1(:,:))) ! G(i w) + G(-i w)
       ! [0.5(b-a)+a,b]
       weval=interval_tmp_coord2(2*iinterval  ,ifreq)
       wweight_eval=interval_tmp_weight2(2*iinterval  ,ifreq)
       weval_cpx=weval*im
       call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_tmp1)
       G_tmp1(:,:)=G_tmp1(:,:)*Exp(im*tcoord(itau)*weval)
       G_test(:,:)=G_test(:,:)+im*wweight_eval*(G_tmp1(:,:)+conjg(G_tmp1(:,:))) ! G(i w) + G(-i w)
      enddo
      G_test=G_test/(2d0*pi)
      G_test(:,:)=G_test(:,:)-interval_vals(itau,:,:)
      do ibas=1,nBas
       do jbas=1,nBas
        err_G_set=err_G_set+abs(G_test(ibas,jbas))
       enddo
      enddo
     enddo
     ! If we entered here, we always retain the new 2 halves of the original interval
     ! If diff. is small we keep the small intervals, but deactivate these two intervals
     err_G_set=err_G_set/(nBas*nBas*ngrid)
     if(err_G_set<thrs_interval) then
      interval_todo2(2*iinterval-1)=.false.
      interval_todo2(2*iinterval  )=.false.
     endif
    else ! The interval is inactive (we store it at the odd row and set to 0 the even row)
     interval_tmp_weight2(2*iinterval-1,:)=interval_tmp_weight(iinterval,:)
     interval_tmp_coord2(2*iinterval-1,:)=interval_tmp_coord(iinterval,:)
     interval_r2(2*iinterval-1,:)=interval_r(iinterval,:)
     interval_todo2(2*iinterval-1)=interval_todo(iinterval)
     interval_tmp_weight2(2*iinterval ,:)=0d0
     interval_tmp_coord2(2*iinterval  ,:)=0d0
     interval_r2(2*iinterval ,:)=0d0
     interval_todo2(2*iinterval)=.false.
    endif 
   enddo
   ! Find the intervals that survived
   nintervals=0
   do iinterval=1,nintervals_twice
    if( abs(interval_r2(iinterval,1))<1d-8 .and. abs(interval_r2(iinterval,2))<1d-8 ) then
     cycle
    else
     nintervals=nintervals+1 
    endif
   enddo
   ! Save in interval_tmp_weight, interval_tmp_coord, interval_r, and interval_todo (i.e., the intervals that survived)
   deallocate(interval_tmp_weight) 
   deallocate(interval_tmp_coord) 
   deallocate(interval_r)
   deallocate(interval_todo)
   allocate(interval_tmp_weight(nintervals,n01))
   allocate(interval_tmp_coord(nintervals,n01))
   allocate(interval_r(nintervals,2))
   allocate(interval_todo(nintervals))
   jinterval=1
   do iinterval=1,nintervals_twice
    if( abs(interval_r2(iinterval,1))<1d-8 .and. abs(interval_r2(iinterval,2))<1d-8 ) then
     cycle
    else
     interval_tmp_weight(jinterval,:)=interval_tmp_weight2(iinterval,:) 
     interval_tmp_coord(jinterval,:)=interval_tmp_coord2(iinterval,:) 
     interval_r(jinterval,:)=interval_r2(iinterval,:) 
     interval_todo(jinterval)=interval_todo2(iinterval)
     jinterval=jinterval+1 
    endif
   enddo
   deallocate(interval_tmp_weight2) 
   deallocate(interval_tmp_coord2) 
   deallocate(interval_r2)
   deallocate(interval_todo2)
  endif

  ! Transfer from the matrices to the vectors wweight and wcoord to compute the full
  !    G_int(itau_k) = i/(2 pi) int _0 ^Infty [ G(i w) Exp(i tau_k w)  + G(-i w) Exp(-i tau_k w) ] d w
  nfreqs=n01*nintervals
  deallocate(wweight,wcoord)
  allocate(wweight(nfreqs),wcoord(nfreqs))
  jgrid=1
  do iinterval=1,nintervals
   do igrid=1,n01
    wweight(jgrid) = interval_tmp_weight(iinterval,igrid)
    wcoord(jgrid)  =  interval_tmp_coord(iinterval,igrid)
    jgrid=jgrid+1
   enddo
  enddo

  ! Check the global error for G(i w) -> G(i tau)
  err_G_set=0d0
  do itau=1,ngrid
   G_test=czero
   do ifreq=1,nfreqs
    weval_cpx=wcoord(ifreq)*im
    call G_AO_RHF(nBas,nOrb,nO,0d0,cHF,eHF,weval_cpx,G_tmp1)
    G_tmp1(:,:)=G_tmp1(:,:)*Exp(im*tcoord(itau)*wcoord(ifreq))
    G_test(:,:)=G_test(:,:)+im*wweight(ifreq)*(G_tmp1(:,:)+conjg(G_tmp1(:,:))) ! G(i w) + G(-i w)
   enddo
   G_test=G_test/(2d0*pi)
   if(verbose/=0) then ! Print for comparison and debug
    write(*,'(a,f25.15)') ' tcoord_test',tcoord(itau)
    write(*,*) 'Integrated'
    do ibas=1,nbas
     write(*,'(*(f15.8))') G_test(ibas,:)
    enddo
    write(*,*) 'Analytic'
    do ibas=1,nbas
     write(*,'(*(f15.8))') G_set_test(itau,ibas,:)
    enddo
    write(*,*)
   endif
   G_test(:,:)=G_test(:,:)-G_set_test(itau,:,:)
   do ibas=1,nBas
    do jbas=1,nBas
     err_G_set=err_G_set+abs(G_test(ibas,jbas))
    enddo
   enddo
  enddo
  err_G_set=err_G_set/(nBas*nBas*ngrid) ! MAE
  write(*,'(a,f15.8,a,i5,a,i8,a)') ' MAE in reproducing G(itau) [Integrated - Analytic]',err_G_set,' at iter ',iter,' with ',nfreqs,' grid'

  if(iter>0) then
   if(abs(err_G_set_old-err_G_set)<1d-8 .or. nintervals_old==nintervals) then
    write(*,'(a,e15.8,a,e15.8)') ' Reducing the MAE accepted for each interval integrations form ',thrs_interval,' to ',1d-1*thrs_interval 
    thrs_interval=1d-1*thrs_interval
    interval_todo(:)=.true.
   endif
  endif
  err_G_set_old=err_G_set
  nintervals_old=nintervals

  all_deactivated=.true.
  do iinterval=1,nintervals
   if(interval_todo(iinterval)) all_deactivated=.false. 
  enddo
  if(all_deactivated) then
   write(*,'(a,e15.8,a,e15.8)') ' Reducing the MAE accepted for each interval integrations form ',thrs_interval,' to ',1d-1*thrs_interval 
   thrs_interval=1d-1*thrs_interval
   interval_todo(:)=.true.
  endif

  if(err_G_set<max_err_G_set) exit

  if(iter==niter_max) exit

 enddo
 write(*,*) 

 deallocate(interval_vals)
 deallocate(interval_tmp_weight)
 deallocate(interval_tmp_coord)
 deallocate(interval_r,interval_todo)
 deallocate(tcoord)

!deallocate(wweight,wcoord)
!allocate(weight_01(1),coord_01(1))
!nfreqs=800
!kind_int = 1
!lim_inf = 0d0; lim_sup = 1d0;
!alpha = 0d0;  beta  = 0d0;
!allocate(wweight(nfreqs),wcoord(nfreqs))
!call cgqf(nfreqs,kind_int,alpha,beta,lim_inf,lim_sup,wcoord,wweight)
!wweight(:)=wweight(:)/((1d0-wcoord(:))**2d0)
!wcoord(:)=wcoord(:)/(1d0-wcoord(:))

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
 deallocate(weight_01,coord_01)
 deallocate(wweight,wcoord)

 ! Recover eHF initial values and deallocate arrays

 eHF(:) = eHF(:)+chem_pot

 deallocate(G_test,G_tmp1)
 deallocate(G_set_test)

end subroutine

