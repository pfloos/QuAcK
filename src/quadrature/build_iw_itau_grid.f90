
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

  logical                        :: file_exists
  logical                        :: all_deactivated
  logical,allocatable            :: interval_todo(:)
  logical,allocatable            :: interval_todo2(:)

  integer                        :: iunit=311
  integer                        :: iter
  integer                        :: itau
  integer                        :: ifreq
  integer                        :: ibas,jbas
  integer                        :: igrid,jgrid
  integer                        :: iinterval,jinterval
  integer                        :: kind_int
  integer                        :: ntimes_max,nfreqs_max
  integer                        :: n01,ngrid,nintervals,nintervals_old,nintervals_twice,niter_max

  double precision               :: error
  double precision               :: ran_num
  double precision               :: thrs_tnorm
  double precision               :: thrs_interval
  double precision               :: chem_pot
  double precision               :: max_err_G_set,err_G_set,err_G_set_old
  double precision               :: teval,tstep,teval_plus,teval_minus,teval_2plus,teval_2minus,grad_teval
  double precision               :: max_weval
  double precision               :: max_teval,min_teval,tweight_eval
  double precision               :: alpha,beta,lim_inf,lim_sup

  double precision,allocatable   :: coord_01(:),weight_01(:)
  double precision,allocatable   :: wweight(:),wcoord(:)
  double precision,allocatable   :: tweight(:),tcoord(:)
  double precision,allocatable   :: interval_r(:,:)
  double precision,allocatable   :: interval_r2(:,:)
  double precision,allocatable   :: interval_tmp_weight(:,:),interval_tmp_coord(:,:)
  double precision,allocatable   :: interval_tmp_weight2(:,:),interval_tmp_coord2(:,:)

  complex*16                     :: G_4vals(4)
  complex*16,allocatable         :: G_test(:,:),G_tmp1(:,:)
  complex*16,allocatable         :: interval_vals(:,:)

! Output variables

  integer,intent(inout)          :: ntimes
  integer,intent(inout)          :: nfreqs
  double precision,intent(inout) :: eHF(nOrb)

!------------------------------------------------------------------------
! Build iw and itau grids 
!------------------------------------------------------------------------

 ! Set global variables
 ntimes_max=2000
 nfreqs_max=2000
 ngrid=40
 thrs_tnorm=1d-5
 min_teval=2d-4  ! Start tau grid at least at 2d-4 because when tau->0 G(i tau) is not continuos [ Lim tau->0+ /= Lim tau->0- ]
 allocate(G_test(nBas,nBas),G_tmp1(nbas,nBas))
 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 eHF(:) = eHF(:)-chem_pot
 inquire(file='max_ntimes_nfreqs', exist=file_exists)
 if(file_exists) then
  open(unit=iunit, form='formatted', file='max_ntimes_nfreqs', status='old')
  read(iunit,*) ntimes_max,nfreqs_max
  close(iunit)
 endif

 ! Set grid from 0 to 1
 n01 = 100 ! From 0 to 1 we always take 100 points
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 allocate(weight_01(n01),coord_01(n01))
 call cgqf(n01,kind_int,alpha,beta,lim_inf,lim_sup,coord_01,weight_01)

 ! Use smallest exponent to find max_teval
 tstep=1d0
 teval=1d0
 error=1d0
 iter=0
 do 
  if(abs(Exp(eHF(nO)*teval)-thrs_tnorm)<error) then
   tstep=1d-1*tstep
   error=1d-1*error
   teval=teval+tstep
  else
   teval=teval+tstep
  endif
  iter=iter+1
  if(iter>10000) exit
  if(abs(Exp(eHF(nO)*teval)-thrs_tnorm)<1d-2) exit
 enddo
 iter=0
 do
  teval_2plus =abs(Exp(eHF(nO)*(teval+2d-3))-thrs_tnorm)
  teval_plus  =abs(Exp(eHF(nO)*(teval+1d-3))-thrs_tnorm)
  teval_minus =abs(Exp(eHF(nO)*(teval-1d-3))-thrs_tnorm)
  teval_2minus=abs(Exp(eHF(nO)*(teval-2d-3))-thrs_tnorm)
  error       =abs(Exp(eHF(nO)*teval)-thrs_tnorm)
  grad_teval = (-teval_2plus+8d0*teval_plus-8d0*teval_minus+teval_2minus)/(12d0*1d-3)
  teval=teval-error/(grad_teval+1d-10)
  iter=iter+1
  if(iter>1000) exit
  if(error<1d-8) exit
 enddo
 max_teval=teval
 max_weval=5d6
 write(*,'(a,2f20.10)') '  Maximum t value and error : ',max_teval,error
 write(*,'(a, f20.10)') '  Minimum t value           : ',min_teval
 write(*,'(a, f20.10)') '  Maximum w value           : ',max_weval
 write(*,*)

 ! Find tweight and tcoord with an adaptative quadrature
 iter=-1
 nintervals=1
 thrs_interval=1d-3
 niter_max=40
 max_err_G_set=1d-3
 ! Build the 40 trial freqs from 0 to max_weval
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
 allocate(tweight(1),tcoord(1))
 allocate(interval_vals(ngrid,2))
 allocate(interval_tmp_weight(nintervals,n01))
 allocate(interval_tmp_coord(nintervals,n01))
 allocate(interval_r(nintervals,2))
 allocate(interval_todo(nintervals))
 iinterval=1
 interval_todo(:)=.true.
 interval_r(iinterval,1)=min_teval ! Start tau grid at least at 2d-4 because when tau->0 G(i tau) is not continuos [ Lim tau->0+ /= Lim tau->0- ]
 interval_r(iinterval,2)=max_teval ! From the scan, when Exp[-eHF_HOMO t] = 1d-5
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
   do iinterval=1,nintervals ! ith interval
    if(interval_todo(iinterval)) then
     ! For the ith active interval [a,b] compute all G_interval(iw_k) = -i int _a ^b [ G(i tau) Exp(-i tau w_k)  + G(-i tau) Exp(i tau w_k) ] d tau
     interval_vals=czero
     do ifreq=1,ngrid
      G_4vals=czero
      do itau=1,n01
       teval=interval_tmp_coord(iinterval,itau)
       tweight_eval=interval_tmp_weight(iinterval,itau)
       G_4vals(1)=G_4vals(1)-im*tweight_eval*im*Exp(teval*eHF(nO))           *Exp(-im*teval*wcoord(ifreq)) !  Go(i tau)
       G_4vals(2)=G_4vals(2)-im*tweight_eval*im*Exp(teval*(eHF(1)-eHF(nOrb)))*Exp(-im*teval*wcoord(ifreq)) ! ~Xo(i tau)
      enddo
      interval_vals(ifreq,1)=G_4vals(1)
      interval_vals(ifreq,2)=G_4vals(2)
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
     ! Check the difference in G_interval(iw_k) computed as [a,b] vs [a,0.5(b-a)+a] + [0.5(b-a)+a,b]
     err_G_set=0d0
     do ifreq=1,ngrid
      G_4vals=czero
      do itau=1,n01
       ! [a,0.5(b-a)+a]
       teval=interval_tmp_coord2(2*iinterval-1,itau)
       tweight_eval=interval_tmp_weight2(2*iinterval-1,itau)
       G_4vals(1)=G_4vals(1)-im*tweight_eval*im*Exp(teval*eHF(nO))           *Exp(-im*teval*wcoord(ifreq)) !  Go(i tau)
       G_4vals(2)=G_4vals(2)-im*tweight_eval*im*Exp(teval*(eHF(1)-eHF(nOrb)))*Exp(-im*teval*wcoord(ifreq)) ! ~Xo(i tau)
       ! [0.5(b-a)+a,b]
       teval=interval_tmp_coord2(2*iinterval  ,itau)
       tweight_eval=interval_tmp_weight2(2*iinterval  ,itau)
       G_4vals(1)=G_4vals(1)-im*tweight_eval*im*Exp(teval*eHF(nO))           *Exp(-im*teval*wcoord(ifreq)) !  Go(i tau)
       G_4vals(2)=G_4vals(2)-im*tweight_eval*im*Exp(teval*(eHF(1)-eHF(nOrb)))*Exp(-im*teval*wcoord(ifreq)) ! ~Xo(i tau)
      enddo
      G_4vals(:)=G_4vals(:)-interval_vals(ifreq,:)
      err_G_set=err_G_set+abs(G_4vals(1))+abs(G_4vals(2))
     enddo
     ! If we entered here, we always retain the new 2 halves of the original interval
     ! If diff. is small we keep the small intervals, but deactivate these two intervals
     err_G_set=err_G_set/(2*ngrid)
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
  ! Transfer from the matrices to the vectors tweight and tcoord to compute the full Fourier transform (i.e., from -Infty to Infty)
  !    G(iw_k) = -i int _0 ^Infty [ G(i tau) Exp(-i tau w_k)  + G(-i tau) Exp(i tau w_k) ] d tau
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
  ! Check the global error for G(i tau) -> G(i w) and ~Xo(i tau) -> ~Xo(i w)
  err_G_set=0d0
  do ifreq=1,ngrid
   G_4vals=czero
   do itau=1,ntimes                        ! This is Go or ~Xo 
    G_4vals(1)=G_4vals(1)-im*tweight(itau)*im*Exp(tcoord(itau)*eHF(nO))           *Exp(-im*tcoord(itau)*wcoord(ifreq)) !  Go(i tau)
    G_4vals(2)=G_4vals(2)-im*tweight(itau)*im*Exp(tcoord(itau)*(eHF(1)-eHF(nOrb)))*Exp(-im*tcoord(itau)*wcoord(ifreq)) ! ~Xo(i tau)
   enddo
   G_4vals(1)=G_4vals(1)-1d0/(im*wcoord(ifreq)-eHF(nO))
   G_4vals(2)=G_4vals(2)-1d0/(im*wcoord(ifreq)-(eHF(1)-eHF(nOrb)))
   err_G_set=err_G_set+abs(G_4vals(1))+abs(G_4vals(2))
  enddo
  err_G_set=err_G_set/(2*ngrid) ! MAE
  write(*,'(a,f8.4,a,i5,a,i8,a)') ' MAE in reproducing G(iw) and X(iw) [Integrated - Analytic]',err_G_set,' at iter ',iter,' with ',ntimes,' grid'
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
  if(ntimes>ntimes_max) exit
 enddo
 write(*,*)
 deallocate(interval_vals)
 deallocate(interval_tmp_weight)
 deallocate(interval_tmp_coord)
 deallocate(interval_r,interval_todo)

 ! Print tweight and tcoord
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
 deallocate(tcoord,tweight)
 deallocate(wcoord)









 nfreqs = 800 ! From 0 to 1 we always take 100 points
 kind_int = 1
 lim_inf = 0d0; lim_sup = 1d0;
 alpha = 0d0;  beta  = 0d0;
 allocate(wweight(nfreqs),wcoord(nfreqs))
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
 deallocate(wcoord,wweight)


 ! Recover eHF initial values and deallocate arrays
 eHF(:) = eHF(:)+chem_pot

 ! Deallocate arrays 
 deallocate(weight_01,coord_01)
 deallocate(G_test,G_tmp1)

end subroutine
