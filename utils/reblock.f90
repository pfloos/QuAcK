!--------------------------------------------------------------------!
! REBLOCK                                                            !
! =======                                                            !
! Neil Drummond, 8.2005                                              !
! (Based on earlier MDT reblock utility for old-style .hist files)   !
!                                                                    !
! This utility performs a statistical analysis of the raw QMC data   !
! held in the .hist file.  The reblocking procedure is               !
! necessary in order to obtain reliable statistical error bars for   !
! the mean values of serially correlated data.  Please refer to the  !
! CASINO manual for further information about the procedure.         !
!                                                                    !
! Changes                                                            !
! NDD  9.05  Data files renamed vmc.hist, etc.  Check for old format.!
!            Other new checks on files.                              !
! NDD  9.05  Changed format of .hist files (again).                  !
! AB   5.06  Add reblocking analysis for future walking estimators.  !
! AB  11.07  Add reblocking analysis of forces and introduce a line  !
!            break in the .hist file to allow storing many items.    !
! NDD 05.08  Rearranged output, to put important stuff at end.       !
! NDD 05.10  Allow for FISQ data in qmc.hist.                        !
!--------------------------------------------------------------------!


MODULE stats_calcs
!-------------------------------------------------------------------!
! A collection of subroutines for performing various statistical    !
! analyses of data.                                                 !
!-------------------------------------------------------------------!
 IMPLICIT NONE


CONTAINS


 SUBROUTINE compute_stats_unweighted(want_skew_kurt,n,data_arr,av,var,skew, &
  &kurt,max_val,min_val)
!-------------------------------------------------------------------!
! Compute mean, variance, skewness, kurtosis and max and min of a   !
! set of data.                                                      !
!-------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  DOUBLE PRECISION,INTENT(in) :: data_arr(n)
  LOGICAL,INTENT(in) :: want_skew_kurt
  DOUBLE PRECISION,INTENT(out) :: av,var,skew,kurt,max_val,min_val
  INTEGER i
  DOUBLE PRECISION sum_delta_x2,sum_delta_x3,sum_delta_x4

  if(n<2)then
   write(6,*)'Can''t compute variance with fewer than two points.'
   stop
  endif

! Compute average.
  av=sum(data_arr(1:n))/dble(n)

! Compute max and min.
  max_val=maxval(data_arr(1:n))
  min_val=minval(data_arr(1:n))

  if(want_skew_kurt)then
! Compute variance, skewness and kurtosis.
   sum_delta_x2=0.d0  ;  sum_delta_x3=0.d0  ;  sum_delta_x4=0.d0
   do i=1,n
    sum_delta_x2=sum_delta_x2+(data_arr(i)-av)**2
    sum_delta_x3=sum_delta_x3+(data_arr(i)-av)**3
    sum_delta_x4=sum_delta_x4+(data_arr(i)-av)**4
   enddo ! i
   var=sum_delta_x2/dble(n-1)
   if(var>0.d0)then
    skew=((sqrt(dble(n-1))*dble(n))/dble(n-2))*sum_delta_x3/sum_delta_x2**1.5d0
    kurt=((dble(n+1)*dble(n)*dble(n-1))/(dble(n-2)*dble(n-3)))*&
     &sum_delta_x4/sum_delta_x2**2-&
     &((dble(n-1)*dble(n-1))/(dble(n-2)*dble(n-3)))*3.d0
   else
    skew=0.d0
    kurt=0.d0
   endif
  else
! Compute variance.
   sum_delta_x2=0.d0
   do i=1,n
    sum_delta_x2=sum_delta_x2+(data_arr(i)-av)**2
   enddo ! i
   var=sum_delta_x2/dble(n-1)
   skew=0.d0
   kurt=0.d0
  endif ! want_skew_kurt

 END SUBROUTINE compute_stats_unweighted


 SUBROUTINE reblock_unweighted(no_pts,data_array,block_length,av, &
  &std_err,delta_std_err)
!---------------------------------------------------------------!
! Compute the unweighted average of the data, and calculate the !
! error bar for a given block length.                           !
!---------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: no_pts,block_length
  DOUBLE PRECISION,INTENT(in) :: data_array(:)
  DOUBLE PRECISION,INTENT(out) :: av,std_err,delta_std_err
  INTEGER i,k,no_blocks,no_pts_in_last_block,j
  DOUBLE PRECISION last_block_weight,var,tot_weight, &
   &tot_weight_sq,block_av,red_tot_weight,rec_block_length

! Compute average of data.
  av=sum(data_array(1:no_pts))/dble(no_pts)

! Number of blocks.
  no_blocks=no_pts/block_length
  rec_block_length=1.d0/dble(block_length)

! Evaluate the sum of the squares of the deviations from the average.
! Weight the last, incomplete block by its size as a fraction of the others.
  var=0.d0
  k=0
  do i=1,no_blocks
   block_av=0.d0
   do j=1,block_length
    k=k+1
    block_av=block_av+data_array(k)
   enddo ! j
   block_av=block_av*rec_block_length
   var=var+(block_av-av)**2
  enddo ! i
  block_av=0.d0
  no_pts_in_last_block=0
  do
   k=k+1
   if(k>no_pts)exit
   no_pts_in_last_block=no_pts_in_last_block+1
   block_av=block_av+data_array(k)
  enddo ! k
  last_block_weight=dble(no_pts_in_last_block)*rec_block_length
  if(no_pts_in_last_block>0)then
   block_av=block_av/dble(no_pts_in_last_block)
   var=var+(block_av-av)**2*last_block_weight
  endif ! last block nonzero

! Evaluate variance, standard error in mean and error in standard error.
  tot_weight=dble(no_blocks)+last_block_weight
  tot_weight_sq=dble(no_blocks)+last_block_weight**2
  red_tot_weight=tot_weight-tot_weight_sq/tot_weight
  var=var/red_tot_weight
  std_err=sqrt(var/tot_weight)
  if(tot_weight>1.d0)then
   delta_std_err=std_err/sqrt(2.d0*(tot_weight-1.d0))
  else
   delta_std_err=0.d0
  endif

 END SUBROUTINE reblock_unweighted


 SUBROUTINE reblock_weighted(no_pts,data_array,weight_array,block_length, &
  &av,std_err,delta_std_err)
!--------------------------------------------------------------!
! Compute the weighted average of the data, and calculate the  !
! error bar for a given block length.                          !
!--------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: no_pts,block_length
  DOUBLE PRECISION,INTENT(in) :: weight_array(no_pts),data_array(no_pts)
  DOUBLE PRECISION,INTENT(out) :: av,std_err,delta_std_err
  INTEGER i,k,no_blocks,no_pts_in_last_block,j
  DOUBLE PRECISION var,tot_weight, &
   &tot_weight_sq,block_av,red_tot_weight,block_weight, &
   &eff_no_blocks

! Compute average of data.
  av=0.d0
  tot_weight=0.d0
  do i=1,no_pts
   av=av+data_array(i)*weight_array(i)
   tot_weight=tot_weight+weight_array(i)
  enddo ! i
  av=av/tot_weight

! Number of blocks
  no_blocks=no_pts/block_length

! Evaluate the sum of the squares of the deviations from the average.
! Last, incomplete block has fewer data points and hence a smaller weight.
  var=0.d0
  tot_weight_sq=0.d0
  k=0
  do i=1,no_blocks
   block_av=0.d0
   block_weight=0.d0
   do j=1,block_length
    k=k+1
    block_av=block_av+data_array(k)*weight_array(k)
    block_weight=block_weight+weight_array(k)
   enddo ! j
   block_av=block_av/block_weight
   var=var+(block_av-av)**2*block_weight
   tot_weight_sq=tot_weight_sq+block_weight**2
  enddo ! i
  block_av=0.d0
  block_weight=0.d0
  no_pts_in_last_block=0
  do
   k=k+1
   if(k>no_pts)exit
   no_pts_in_last_block=no_pts_in_last_block+1
   block_av=block_av+data_array(k)*weight_array(k)
   block_weight=block_weight+weight_array(k)
  enddo ! k
  if(no_pts_in_last_block>0)then
   block_av=block_av/block_weight
   var=var+(block_av-av)**2*block_weight
   tot_weight_sq=tot_weight_sq+block_weight**2
  endif ! last block nonzero

! Evaluate variance, standard error in mean and error in standard error.
  red_tot_weight=tot_weight-tot_weight_sq/tot_weight
  var=var/red_tot_weight

  eff_no_blocks=dble(no_blocks)+dble(no_pts_in_last_block)/dble(block_length)

  std_err=sqrt(var/eff_no_blocks)
  if(eff_no_blocks>1.d0)then
   delta_std_err=std_err/sqrt(2.d0*(eff_no_blocks-1.d0))
  else
   delta_std_err=0.d0
  endif

 END SUBROUTINE reblock_weighted


 SUBROUTINE correlation_time(n,Odata,Otau,Otau_err,Oave_in,Ovar_in)
!------------------------------------------------------------------------!
! Obtain correlation time from a set of data                             !
!------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: n
 DOUBLE PRECISION,INTENT(in) :: Odata(n)
 DOUBLE PRECISION,INTENT(in),OPTIONAL :: Oave_in,Ovar_in
 DOUBLE PRECISION,INTENT(out) :: Otau,Otau_err
 DOUBLE PRECISION Oave,Oave2,O2ave,Ovar,invOvar,Oacorr,ri,invn
 DOUBLE PRECISION,PARAMETER :: tol=1.d-100
 INTEGER i,sqrtn

 Otau=-1.d0 ; Otau_err=-1.d0
 if(n<10)return
 invn=1.d0/dble(n)
 sqrtn=nint(sqrt(dble(n)))

! <O>, <O>**2, <O**2>, variance
 if(present(Oave_in))then
  Oave=Oave_in
 else
  Oave=sum(Odata)*invn
 endif
 if(present(Ovar_in))then
  Ovar=Ovar_in*invn*(n-1)
 else
  Oave2=Oave**2
  O2ave=sum(Odata**2)*invn
  Ovar=O2ave-Oave2
 endif
 if(Ovar<tol)return
 invOvar=1.d0/Ovar

! Autocorrelation for i<=cut-off -> tau
 Otau=1.d0
 do i=1,n-1
  Oacorr=sum((Odata(1:n-i)-Oave)*(Odata(1+i:n)-Oave))*invOvar/dble(n-i)
  Otau=Otau+2*Oacorr
  if(i>=nint(3*Otau))then
   ri=dble(i) ; exit
  endif
 enddo

! Error in tau
 Otau_err=Otau*sqrt((4*ri+2.d0)*invn)

 END SUBROUTINE correlation_time


END MODULE stats_calcs


MODULE analysis
!-------------------------------------------------------------!
! Miscellaneous subroutines for reading & analysing the data. !
!-------------------------------------------------------------!
 USE stats_calcs
 IMPLICIT NONE

! Tags for the columns of the data file, specifying where each data item
! is held.  If a tag is negative, the data item isn't present.
 INTEGER tag_step,tag_energy,tag_etotalt,tag_esqr,tag_popavgsqr,tag_K,tag_T, &
  &tag_fisq,tag_Ewald,tag_local,tag_nonlocal,tag_short,tag_long,tag_cppei, &
  &tag_cppe,tag_cppee,tag_masspol,tag_massvel,tag_darwinen,tag_darwinee, &
  &tag_retard,tag_weight,tag_nconf,tag_eref,tag_ebest,tag_acc,tag_teff, &
  &tag_dipole1,tag_dipole2,tag_dipole3,tag_dipole_sq,tag_contact_den, &
  &tag_future0,tag_future1,tag_future2,tag_future3,tag_future4,tag_future5, &
  &tag_future6,tag_future7,tag_future8,tag_future9,tag_future10

! Number of columns of data in .hist file.
 INTEGER no_cols_qmc

! Title of .hist file
 CHARACTER(72) title

! File version number
 INTEGER version

! CASINO input keywords: interaction type and basis type.
 CHARACTER(20) interaction,atom_basis_type

! Do we have Ewald/Coulomb interaction? Do we have MPC interaction?
 LOGICAL coul_mpc,coul_ewald

! Ion-ion energy
 DOUBLE PRECISION constant_energy

! Total number of electrons; no. of atoms per prim cell; no. primitive cells
 INTEGER netot,nbasis,npcells

! Number of parts in simulation cell:
!  =npcells for periodic systems
!  =netot for electron gas
!  =1 otherwise
 INTEGER nparts_per_simcell

! Is the system periodic?
 LOGICAL isperiodic

! QMC method used
 CHARACTER(3) qmc_method

! Number of lines of data.
 INTEGER Nlines

! Number of equilibration lines.
 INTEGER Nequil

! Name of .hist file.
 CHARACTER(8) filename

! Number of initial lines to discard.
 INTEGER Nskip

! Array with the hist data from the files.
 DOUBLE PRECISION,ALLOCATABLE :: data_array(:,:)

! Energy units
 CHARACTER(15) e_units

! Units conversion: a.u.->eV and a.u.->kcal.
 DOUBLE PRECISION,PARAMETER :: htoev=27.2113962d0,htokcal=627.507541278d0

! Are forces to be calculated?
 INTEGER iion,iaxis,item,nitot_forces,naxis_forces,nitem_forces,&
  &nitot_max_forces
 INTEGER,ALLOCATABLE :: tag_forces(:,:,:)
 DOUBLE PRECISION,ALLOCATABLE :: forces_array(:,:)
 LOGICAL forces

 ! Use weights when calculating average energy, etc.
 LOGICAL,PARAMETER :: use_weights=.true.


CONTAINS


 SUBROUTINE read_header(io,dmc)
!----------------------------------------------------------------------!
! Read in the data in the .hist file header and count the lines, etc.  !
!----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: io
  LOGICAL,INTENT(in) :: dmc
  INTEGER ierr,i,s,isper_flag,nbreak,ialloc
  CHARACTER(1) temp
  CHARACTER(72) datastring
  CHARACTER(500) checkline
  LOGICAL,PARAMETER :: verbose=.false.

! Check we don't have an old-style .hist file.
  rewind(io)
  read(io,'(a)',iostat=ierr)checkline
  call check_ierr(ierr)
  checkline=adjustl(checkline)
  if(index(checkline,'Block')>0)then
   write(6,*)'You appear to be analyzing a CASINO version 1 vmc.hist file. &
    &Please use the'
   write(6,*)'UPDATE_HIST utility to update it to the new format.'
   stop
  endif ! Old-style vmc.hist
  if(index(checkline,'#')==0)then
   write(6,*)'Your data file does not seem to start with a header.  &
    &This may be because you'
   write(6,*)'are using an old-format file.  If this is the case then &
    &please use UPDATE_HIST'
   write(6,*)'to update your file.'
   stop
  endif ! No header.

! Count the data lines.  Ignore comments.
  rewind(io)
  Nlines=0
  Nequil=0
  forces=.false.
  do
   read(io,'(a)',iostat=ierr)datastring
   if(ierr>0)then
    write(6,*)'Error reading data file.'
    stop
   endif
   if(ierr<0)exit
   if(index(datastring,'#')==0)then
    Nlines=Nlines+1
   else
    if(trim(adjustl(datastring))=='#### START STATS')Nequil=Nlines
   endif ! Line not a comment.
   if(.not.forces)then
    if(index(datastring,'FOR')>0)forces=.true. ! atomic forces present
   endif
  enddo ! lines
  if(dmc.and.Nequil==0)Nequil=Nlines
  rewind(io)
  if(verbose)then
   if(Nlines/=1)then
    write(6,*)'There are '//trim(i2s(Nlines))//' lines of data in ' &
     &//trim(filename)//'.'
   else
    write(6,*)'There is 1 line of data in '//trim(filename)//'.'
   endif ! Singular / plural
   if(Nequil>1)then
    write(6,*)'Of these, '//trim(i2s(Nlines))//' lines are marked as &
     &equilibration data.'
   elseif(Nequil==1)then
    write(6,*)'Of these, 1 line is marked as equilibration data.'
   else
    write(6,*)'No data are marked as equilibration data.'
   endif ! Nequil
  endif ! verbose
  if(Nlines<2)then
   write(6,*)'There are less than two lines of data in '//trim(filename)//'.'
   write(6,*)'One cannot obtain error bars with fewer than 2 data points.'
   stop
  endif

! Get title.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,'(a)',iostat=ierr)datastring
  call check_ierr(ierr)
  s=index(datastring,'#')
  if(s>0)then
   title=datastring(s+1:len_trim(datastring))
  else
   write(6,*)'Header line does not have a "#" in front.  Stopping.'
   stop
  endif
  title=adjustl(title)
  if(verbose)write(6,*)'Title: '//trim(title)

! Get version number.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,version
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)write(6,*)'File version number is '//trim(i2s(version))//'.'
  if(version/=1)then
   write(6,*)'Version number of '//trim(filename)//' must be 1.'
   stop
  endif ! version/=1

! Get QMC method.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,qmc_method
  call check_ierr(ierr)
  call check_hash(temp)
  qmc_method=adjustl(qmc_method)
  if(verbose)write(6,*)'The data were generated using '//trim(qmc_method)//'.'
  if(trim(qmc_method)/='VMC'.and.trim(qmc_method)/='DMC')then
   write(6,*)'Method in '//trim(filename)//' should be either VMC or DMC.'
   stop
  endif ! method
  if(trim(filename)=='vmc.hist'.and.trim(qmc_method)/='VMC')then
   write(6,*)'Warning: you appear to have non-VMC data in a file called &
    &vmc.hist.'
   write(6,*)
  endif
  if(trim(filename)=='dmc.hist'.and.trim(qmc_method)/='DMC')then
   write(6,*)'Warning: you appear to have non-DMC data in a file called &
    &dmc.hist.'
   write(6,*)
  endif

! Get interaction-type (interaction).
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,interaction
  call check_ierr(ierr)
  call check_hash(temp)
  coul_ewald=.false. ; coul_mpc=.false.
  select case(trim(interaction))
  case('none','coulomb','ewald','mpc','ewald_mpc','mpc_ewald','manual')
   continue
  case('1') ; interaction='default'
  case('2') ; interaction='mpc'
  case('3') ; interaction='ewald_mpc'
  case('4') ; interaction='mpc_ewald'
  case default
   write(6,*)'Value of INTERACTION=',trim(interaction),' not recognized. &
    &Stopping.'
   stop
  end select
  select case(trim(interaction))
  case('none') ; continue
  case('coulomb','ewald','default','manual') ; coul_ewald=.true.
  case('mpc') ; coul_mpc=.true.
  case('ewald_mpc','mpc_ewald') ; coul_ewald=.true. ; coul_mpc=.true.
  end select
  if(verbose)write(6,*)'The value of the interaction parameter is ',&
   &trim(interaction),'.'

! Get constant (ion-ion) energy.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,constant_energy
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)write(6,*)'Have got constant energy component.'

! Get total number of electrons.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,netot
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)then
   if(netot/=1)then
    write(6,*)'There are '//trim(i2s(netot))//' particles in the simulation.'
   else
    write(6,*)'There is 1 particle in the simulation.'
   endif
  endif ! verbose
  if(netot<1)then
   write(6,*)'Should be more than one particle!'
   stop
  endif

! Get number of atoms per primitive cell.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,nbasis
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)then
   if(nbasis/=1)then
    write(6,*)'The primitive cell contains '//trim(i2s(nbasis))//' atoms.'
   else
    write(6,*)'The primitive cell contains 1 atom.'
   endif
  endif ! verbose
  if(nbasis<0)then
   write(6,*)'There should be at least zero atoms...'
   stop
  endif

! Get number of primitive cells.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,npcells
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)then
   if(npcells/=1)then
    write(6,*)'There are '//trim(i2s(npcells))//' primitive cells.'
   else
    write(6,*)'There is 1 primitive cell.'
   endif
  endif ! verbose
  if(npcells<1)then
   write(6,*)'There should be at least one primitive cell.'
   stop
  endif

! When forces are present, allocate force array.
  if(forces)then
   nitot_max_forces=nbasis*npcells
   allocate(tag_forces(22,3,nitot_max_forces),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)'Force array allocation problem.'
    stop
   endif ! ialloc/=0
  endif ! forces

! Basis-type keyword.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,atom_basis_type
  call check_ierr(ierr)
  call check_hash(temp)
  select case(trim(atom_basis_type))
  case('0') ; atom_basis_type='none'
  case('1') ; atom_basis_type='plane-wave'
  case('2') ; atom_basis_type='gaussian'
  case('3') ; atom_basis_type='numerical'
  case('4') ; atom_basis_type='blip'
  case('5') ; atom_basis_type='non_int_he'
  case default
   continue
  end select
  if(verbose)write(6,*)'The value of the atom_basis_type parameter is ' &
   &//trim(atom_basis_type)//'.'

! Get periodicity.
  read(io,*,iostat=ierr)
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,isper_flag
  call check_ierr(ierr)
  call check_hash(temp)
  if(isper_flag==1)then
   isperiodic=.true. ; if(verbose)write(6,*)'The system is periodic.'
   select case(trim(interaction))
   case('default','coulomb') ; interaction='ewald'
   end select
  elseif(isper_flag==0)then
   isperiodic=.false. ; if(verbose)write(6,*)'The system is not periodic.'
   select case(trim(interaction))
   case('default','ewald') ; interaction='coulomb'
   case('mpc','mpc_ewald','ewald_mpc')
    write(6,*)'Interaction type should be ''coulomb'' or ''none'' for finite &
     &systems. Contradiction in header.'
    stop
   end select
  else
   write(6,*)'Periodicity flag must be 0 or 1.'
   stop
  endif ! periodicity.

! Get number of data columns.  Increase it by 1, since the line-numbers will
! also be read.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  read(io,*,iostat=ierr)temp,no_cols_qmc
  call check_ierr(ierr)
  call check_hash(temp)
  if(verbose)then
   if(no_cols_qmc/=1)then
    write(6,*)'There are '//trim(i2s(no_cols_qmc))//' columns of data in ' &
     &//trim(filename)//'.'
   else
    write(6,*)'There is 1 column of data in '//trim(filename)//'.'
   endif ! Singular/plural
  endif ! verbose
  if(no_cols_qmc<1)then
   write(6,*)'No data to analyse.  Stopping.'
   stop
  endif
  no_cols_qmc=no_cols_qmc+1
! Account for line breaking as the maximum number of items per line is 25
  nbreak=no_cols_qmc/25
  if(modulo(no_cols_qmc,25)>0)nbreak=nbreak+1
  Nlines=Nlines/nbreak
  Nequil=Nequil/nbreak

! Get items in .hist file
  tag_step=1         ! Move number
  tag_energy=-1      ! Total energy
  tag_etotalt=-1     ! Alternative total energy
  tag_esqr=-1        ! Square of total energy
  tag_popavgsqr=-1   ! Square of population average over total energy
  tag_K=-1           ! KEI kinetic-energy estimator
  tag_T=-1           ! TI kinetic-energy estimator
  tag_fisq=-1        ! FISQ kinetic-energy estimator
  tag_Ewald=-1       ! 1/r or Ewald e-e interaction
  tag_local=-1       ! Local electron-ion energy
  tag_nonlocal=-1    ! Nonlocal electron-ion energy
  tag_short=-1       ! Short-range part of MPC
  tag_long=-1        ! Long-range part of MPC
  tag_cppei=-1       ! Electron-ion CPP term
  tag_cppe=-1        ! Electron CPP term
  tag_cppee=-1       ! Electron-electron CPP term
  tag_masspol=-1     ! Mass-polarization term
  tag_future0=-1     ! Future-walking estimator
  tag_future1=-1     ! "
  tag_future2=-1     ! "
  tag_future3=-1     ! "
  tag_future4=-1     ! "
  tag_future5=-1     ! "
  tag_future6=-1     ! "
  tag_future7=-1     ! "
  tag_future8=-1     ! "
  tag_future9=-1     ! "
  tag_future10=-1    ! "
  tag_massvel=-1     ! Mass-velocity term
  tag_darwinen=-1    ! Darwin e-n term
  tag_darwinee=-1    ! Darwin e-e term
  tag_retard=-1      ! Retardation term.
  tag_weight=-1      ! Total weight of configs
  tag_nconf=-1       ! Number of configs
  tag_eref=-1        ! Reference energy
  tag_ebest=-1       ! Best estimate of energy
  tag_acc=-1         ! Acceptance ratio
  tag_teff=-1        ! Effective time step
  tag_dipole1=-1     ! Electric dipole moment
  tag_dipole2=-1     ! "        "      "
  tag_dipole3=-1     ! "        "      "
  tag_dipole_sq=-1   ! "        "      "
  tag_contact_den=-1 ! Electron-positron contact density
  if(forces)then
   tag_forces(1:22,1:3,1:nitot_max_forces)=-1
   nitem_forces=0 ; naxis_forces=0 ; nitot_forces=0
  endif ! forces

  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)
  do i=2,no_cols_qmc
   read(io,*,iostat=ierr)temp,datastring
   call check_ierr(ierr)
   call check_hash(temp)
   datastring=adjustl(datastring)
   if(trim(datastring)=='ETOT')then
    call check_tag_free(tag_energy)
    tag_energy=i
   elseif(trim(datastring)=='ETOTALT')then
    call check_tag_free(tag_etotalt)
    tag_etotalt=i
   elseif(trim(datastring)=='ESQR')then
    call check_tag_free(tag_esqr)
    tag_esqr=i
   elseif(trim(datastring)=='POPAVGSQR')then
    call check_tag_free(tag_popavgsqr)
    tag_popavgsqr=i
   elseif(trim(datastring)=='KEI')then
    call check_tag_free(tag_K)
    tag_K=i
   elseif(trim(datastring)=='TI')then
    call check_tag_free(tag_T)
    tag_T=i
   elseif(trim(datastring)=='FISQ')then
    call check_tag_free(tag_fisq)
    tag_fisq=i
   elseif(trim(datastring)=='EWALD')then
    call check_tag_free(tag_Ewald)
    tag_Ewald=i
   elseif(trim(datastring)=='LOCAL')then
    call check_tag_free(tag_local)
    tag_local=i
   elseif(trim(datastring)=='NONLOCAL')then
    call check_tag_free(tag_nonlocal)
    tag_nonlocal=i
   elseif(trim(datastring)=='SHORT')then
    call check_tag_free(tag_short)
    tag_short=i
   elseif(trim(datastring)=='LONG')then
    call check_tag_free(tag_long)
    tag_long=i
   elseif(trim(datastring)=='CPPEI')then
    call check_tag_free(tag_cppei)
    tag_cppei=i
   elseif(trim(datastring)=='CPPE')then
    call check_tag_free(tag_cppe)
    tag_cppe=i
   elseif(trim(datastring)=='CPPEE')then
    call check_tag_free(tag_cppee)
    tag_cppee=i
   elseif(trim(datastring)=='MASSPOL')then
    call check_tag_free(tag_masspol)
    tag_masspol=i
   elseif(trim(datastring(1:3))=='FOR')then
    call generate_tag_forces(datastring,i)
   elseif(trim(datastring)=='FUTURE0')then
    call check_tag_free(tag_future0)
    tag_future0=i
   elseif(trim(datastring)=='FUTURE1')then
    call check_tag_free(tag_future1)
    tag_future1=i
   elseif(trim(datastring)=='FUTURE2')then
    call check_tag_free(tag_future2)
    tag_future2=i
   elseif(trim(datastring)=='FUTURE3')then
    call check_tag_free(tag_future3)
    tag_future3=i
   elseif(trim(datastring)=='FUTURE4')then
    call check_tag_free(tag_future4)
    tag_future4=i
   elseif(trim(datastring)=='FUTURE5')then
    call check_tag_free(tag_future5)
    tag_future5=i
   elseif(trim(datastring)=='FUTURE6')then
    call check_tag_free(tag_future6)
    tag_future6=i
   elseif(trim(datastring)=='FUTURE7')then
    call check_tag_free(tag_future7)
    tag_future7=i
   elseif(trim(datastring)=='FUTURE8')then
    call check_tag_free(tag_future8)
    tag_future8=i
   elseif(trim(datastring)=='FUTURE9')then
    call check_tag_free(tag_future9)
    tag_future9=i
   elseif(trim(datastring)=='FUTURE10')then
    call check_tag_free(tag_future10)
    tag_future10=i
   elseif(trim(datastring)=='MASSVEL')then
    call check_tag_free(tag_massvel)
    tag_massvel=i
   elseif(trim(datastring)=='DARWINEN')then
    call check_tag_free(tag_darwinen)
    tag_darwinen=i
   elseif(trim(datastring)=='DARWINEE')then
    call check_tag_free(tag_darwinee)
    tag_darwinee=i
   elseif(trim(datastring)=='RETARD')then
    call check_tag_free(tag_retard)
    tag_retard=i
   elseif(trim(datastring)=='WEIGHT')then
    call check_tag_free(tag_weight)
    tag_weight=i
   elseif(trim(datastring)=='NCONF')then
    call check_tag_free(tag_nconf)
    tag_nconf=i
   elseif(trim(datastring)=='EREF')then
    call check_tag_free(tag_eref)
    tag_eref=i
   elseif(trim(datastring)=='EBEST')then
    call check_tag_free(tag_ebest)
    tag_ebest=i
   elseif(trim(datastring)=='ACC')then
    call check_tag_free(tag_acc)
    tag_acc=i
   elseif(trim(datastring)=='TEFF')then
    call check_tag_free(tag_teff)
    tag_teff=i
   elseif(trim(datastring)=='DIPOLE1')then
    call check_tag_free(tag_dipole1)
    tag_dipole1=i
   elseif(trim(datastring)=='DIPOLE2')then
    call check_tag_free(tag_dipole2)
    tag_dipole2=i
   elseif(trim(datastring)=='DIPOLE3')then
    call check_tag_free(tag_dipole3)
    tag_dipole3=i
   elseif(trim(datastring)=='DIPOLESQ')then
    call check_tag_free(tag_dipole_sq)
    tag_dipole_sq=i
   elseif(trim(datastring)=='CONTACT_DEN')then
    call check_tag_free(tag_contact_den)
    tag_contact_den=i
   else
    write(6,*)'Column label not recognised.'
    write(6,*)'Label is: '//trim(datastring)
    stop
   endif ! Label
  enddo ! i
  if(verbose)then
   write(6,*)'Have read in column labels.'
   write(6,*)
  endif ! verbose

! Warn about missing data, etc.
  if(tag_energy<=0)then
   write(6,*)'Warning: total energy data are not present!'
   write(6,*)
  endif
  if(tag_K<=0)then
   write(6,*)'Warning: kinetic energy (K) data are not present!'
   write(6,*)
  endif
  if((tag_short>0.or.tag_long>0).and..not.coul_mpc)then
   write(6,*)'Warning: MPC data are inexplicably present.'
   write(6,*)
  endif
  if(tag_ewald>0.and..not.coul_ewald)then
   write(6,*)'Warning: Ewald data are inexplicably present.'
   write(6,*)
  endif
  if(tag_short>0.and.tag_long<=0)then
   write(6,*)'Warning: only have short-ranged part of MPC interaction.'
   write(6,*)
  endif
  if(tag_short<=0.and.tag_long>0)then
   write(6,*)'Warning: only have long-ranged part of MPC interaction.'
   write(6,*)
  endif

! Read final comment line in header.
  read(io,*,iostat=ierr)temp
  call check_ierr(ierr)
  call check_hash(temp)

 END SUBROUTINE read_header


 SUBROUTINE check_hash(char)
!---------------------------------------------------------------------------!
! This sub is used to check that the 1st char in each header line is a "#". !
!---------------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(1),INTENT(in) :: char
  if(char/='#')then
   write(6,*)'Header line does not have a "#" in front.  Stopping.'
   stop
  endif
 END SUBROUTINE check_hash


 SUBROUTINE check_ierr(ierr,nline)
!------------------------------------------------------!
! Complain if there has been a problem reading a file. !
!------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: ierr
  INTEGER,INTENT(in),OPTIONAL :: nline
  if(ierr/=0)then
   if(present(nline))then
    write(6,*)'Problem reading '//trim(filename)//' at line '//trim(i2s(nline))//'.'
   else
    write(6,*)'Problem reading '//trim(filename)//'.'
   endif
   stop
  endif
 END SUBROUTINE check_ierr


 SUBROUTINE check_tag_free(tag)
!----------------------------------------------!
! Complain if a tag has already been assigned. !
!----------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: tag
  if(tag/=-1)then
   write(6,*)'Tag assigned twice.  Two column labels must be the same.'
   stop
  endif
 END SUBROUTINE check_tag_free


 SUBROUTINE read_data(dmc)
!--------------------------------------------------!
! Read in the raw QMC data from the .hist file.    !
!--------------------------------------------------!
  IMPLICIT NONE
  LOGICAL,INTENT(in) :: dmc
  INTEGER ierr,i,ialloc,nbreak,nleft,in,im
  INTEGER,PARAMETER :: io=8
  CHARACTER(640) char640

! Open the data file.
  open(unit=io,file=trim(filename),status='old',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Sorry, cannot open '//trim(filename)//'.'
   stop
  endif

! Count the columns and rows of data, establish which data are present, etc.
  call read_header(io,dmc)

! Allocate the data array.
  allocate(data_array(Nlines,no_cols_qmc),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation problem (1).'
   stop
  endif

! Read in the data.  Ignore comments.
  i=0
  do
   read(io,'(a)',iostat=ierr)char640
   call check_ierr(ierr,i+1)
   if(index(char640,'#')==0)then
    i=i+1
! When reading from .hist file, account for maximum number of items
! (=25 per line).
    nbreak=no_cols_qmc/25
    nleft=modulo(no_cols_qmc,25)
    im=0
    if(nbreak>0)then
     do in=1,nbreak
      im=in
      read(char640,*,iostat=ierr)data_array(i,(im-1)*25+1:im*25)
      read(io,'(a)',iostat=ierr)char640
     enddo
    endif
    if(nleft>0)then
     read(char640,*,iostat=ierr)data_array(i,im*25+1:no_cols_qmc)
    endif
    call check_ierr(ierr,i+1)
    if(i>=Nlines)exit
   endif ! Line not a comment.
  enddo ! i
  close(io)

 END SUBROUTINE read_data


 SUBROUTINE check_data
!------------------------------------------------------------------------!
! This subroutine checks that the raw data in data_array are consistent. !
! It checks that adding up the energy components gives the total energy, !
! and that the ion-ion energy in the header is correct.  It looks for    !
! Ewald and MPC data and decides which is to be used in the total energy.!
! The number of equilibration steps to be discarded are chosen and the   !
! energy units are selected.                                             !
!------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER i,ierr,units_choice
  DOUBLE PRECISION econst_check,econst_hist,escale,tot_weight
  LOGICAL econst_is_const
! Tolerance for checking that total energy is sum of components.
  DOUBLE PRECISION,PARAMETER :: tol=1.d-6

! Check move numbers
  if(tag_step>0)then
   do i=1,Nlines
    if(nint(data_array(i,tag_step))/=i)then
     write(6,*)'WARNING: iteration number behaves oddly at line ' &
      &//trim(i2s(i))//' in '//trim(filename)//'.'
     write(6,*)
     exit
    endif ! Problem with move number
   enddo ! i
  endif ! tag_step>0

! Check weights.
  if(tag_weight>0)then
   tot_weight=0.d0
   do i=1,Nlines
    if(data_array(i,tag_weight)<0.d0)then
     write(6,*)'Found a negative weight at line '//trim(i2s(i)) &
      &//' of '//trim(filename)//'.'
     stop
    endif ! weight<0
    tot_weight=tot_weight+data_array(i,tag_weight)
   enddo ! i
   if(tot_weight<=0.d0)then
    write(6,*)'Sum of weights is 0.  Stopping.'
    stop
   endif ! total weight=0
   if(.not.use_weights)then
    write(6,*)'Weights are present, but will not be used.'
    write(6,*)
   endif ! weights not to be used.
  endif ! weights present.

! Check that total energy minus KE, e-i pot E, and e-e pot E is a constant:
! the ion-ion energy.
  if(tag_energy>0)then
   econst_is_const=.true.
   do i=1,Nlines
    econst_check=data_array(i,tag_energy)
    if(tag_K>0)econst_check=econst_check-data_array(i,tag_K)
    if(trim(interaction)=='coulomb'.or.trim(interaction)=='ewald'.or.&
     &trim(interaction)=='ewald_mpc'.or.trim(interaction)=='manual')then
     if(tag_ewald>0)econst_check=econst_check-data_array(i,tag_ewald)
    elseif(trim(interaction)=='mpc'.or.trim(interaction)=='mpc_ewald')then
     if(tag_short>0)econst_check=econst_check-data_array(i,tag_short)
     if(tag_long>0)econst_check=econst_check-data_array(i,tag_long)
    endif ! MPC or Ewald present.
    if(tag_local>0)econst_check=econst_check-data_array(i,tag_local)
    if(tag_nonlocal>0)econst_check=econst_check-data_array(i,tag_nonlocal)
    if(tag_cppei>0)econst_check=econst_check-data_array(i,tag_cppei)
    if(tag_cppe>0)econst_check=econst_check-data_array(i,tag_cppe)
    if(tag_cppee>0)econst_check=econst_check-data_array(i,tag_cppee)
    if(i==1)then
     econst_hist=econst_check
    else
     if(abs(econst_hist-econst_check)>tol)then
      write(6,*)'Warning: some component of energy is not accounted for!'
      write(6,*)'First evaluation of ion-ion energy: ',econst_hist
      write(6,*)'Later evaluation of ion-ion energy: ',econst_check
      write(6,*)
      econst_is_const=.false.
      exit
     endif ! constanet_energy not constant.
    endif ! i=1
   enddo ! i
   if(abs(econst_hist)<tol)econst_hist=0.d0
   if(econst_is_const)then
    if(abs(constant_energy-econst_hist)>tol)then
     write(6,*)'Warning: value of ion-ion energy obtained from raw data &
      &differs from the value'
     write(6,*)'in the header.  Missing constant energy component?'
     write(6,*)
    endif ! Difference in constant_energy
   endif ! Can compare constant_energy values
  endif ! Components for check present?

! Check that FISQ=2*TI-KEI.
  if(tag_K>0.and.tag_T>0.and.tag_fisq>0)then
   do i=1,Nlines
    if(abs(data_array(i,tag_fisq)+data_array(i,tag_K) &
     &-2.d0*data_array(i,tag_T))>tol)then
     write(6,*)'Warning: problem with kinetic-energy estimators.  &
      &FISQ /= 2.TI-KEI.'
     write(6,*)' KEI at line '//trim(i2s(i))//' : ',data_array(i,tag_K)
     write(6,*)'  TI at line '//trim(i2s(i))//' : ',data_array(i,tag_T)
     write(6,*)'FISQ at line '//trim(i2s(i))//' : ',data_array(i,tag_fisq)
     write(6,*)
     exit
    endif ! Problem with KE estimators.
   enddo ! i
  endif ! All KE estimators present.

! Find out how many lines are to be skipped.
  if(trim(qmc_method)=='DMC')then
   do
    write(6,*)'There are '//trim(i2s(Nlines))//' lines of data in total.'
    write(6,*)'There are '//trim(i2s(Nequil))//' lines of equilibration data.'
    write(6,*)'How many initial lines of data do you wish to discard?'
    read(5,*,iostat=ierr)Nskip
    if(ierr/=0)Nskip=-1
    if(Nskip<0.or.Nskip>Nlines-2)then
     write(6,*)'Number of lines to skip must be between 0 and '// &
      &trim(i2s(Nlines-2))//'.'
    else
     exit
    endif ! Problem with Nskip
   enddo ! Loop asking for Nskip
   if(Nskip<Nequil)then
    write(6,*)'Warning: equilibration data will be included in the &
     &statistics analysis.'
    write(6,*)
   endif
  else
! Don't skip any lines for VMC.
   Nskip=0
  endif ! QMC method
  write(6,*)trim(i2s(Nlines-Nskip))//' lines of data will be analysed.'
  write(6,*)'These data start on line '//trim(i2s(Nskip+1)) &
   &//' and end on line '//trim(i2s(Nlines))//'.'
  write(6,*)

! Find out how input data is scaled w.r.t. the simulation cell
  if(trim(atom_basis_type)=='none')then
   nparts_per_simcell=netot
  elseif(isperiodic)then
   nparts_per_simcell=npcells
  else
   nparts_per_simcell=1
  endif

! Find out what units are to be used.
  if(trim(atom_basis_type)=='none')then
! Electron gas.
   do
    write(6,*)'Your data are for an electron(-hole) system.'
    write(6,*)'Please select units for your energy data.'
    write(6,*)'Choose one of: (1) au per particle; &
     &(2) eV per particle.'
    read(5,*,iostat=ierr)units_choice
    if(ierr/=0)units_choice=-1
    if(units_choice>=1.and.units_choice<=2)then
     exit
    else
     write(6,*)'Please try again.  Choose a number between 1 and 2.'
     write(6,*)
    endif
   enddo ! choice loop
   if(units_choice==1)then
    escale=1.d0
    e_units='au/particle'
   else
    escale=htoev
    e_units='eV/particle'
   endif ! units choice
  elseif(isperiodic.and.nbasis>0)then
! Periodic system.
   do
    write(6,*)'Your data are for a periodic system with atoms.'
    write(6,*)'Please select units for your energy data.'
    write(6,*)'Choose one of: (1) au per prim cell; (2) eV per &
     &prim cell;'
    write(6,*)'(3) kcal per prim cell; (4) au per atom; (5) eV per atom; &
     &(6) kcal per atom.'
    read(5,*,iostat=ierr)units_choice
    if(ierr/=0)units_choice=-1
    if(units_choice>=1.and.units_choice<=6)then
     exit
    else
     write(6,*)'Please try again.  Choose a number between 1 and 6.'
     write(6,*)
    endif
   enddo ! choice loop
   if(units_choice==1)then
    escale=1.d0
    e_units='au/prim cell'
   elseif(units_choice==2)then
    escale=htoev
    e_units='eV/prim cell'
   elseif(units_choice==3)then
    escale=htokcal
    e_units='kcal/prim cell'
   elseif(units_choice==4)then
    escale=1.d0/dble(nbasis)
    e_units='au/atom'
   elseif(units_choice==5)then
    escale=htoev/dble(nbasis)
    e_units='eV/atom'
   else
    escale=htokcal/dble(nbasis)
    e_units='kcal/atom'
   endif
  elseif(.not.isperiodic.and.nbasis>0)then
! Finite system.
   do
    write(6,*)'Your data are for a finite system with atoms.'
    write(6,*)'Please select units for your energy data.'
    write(6,*)'Choose one of: (1) au; (2) eV; (3) kcal; (4) au per atom; &
     &(5) eV per atom;'
    write(6,*)'(6) kcal per atom.'
    read(5,*,iostat=ierr)units_choice
    if(ierr/=0)units_choice=-1
    if(units_choice>=1.and.units_choice<=6)then
     exit
    else
     write(6,*)'Please try again.  Choose a number between 1 and 6.'
     write(6,*)
    endif
   enddo ! choice loop
   if(units_choice==1)then
    escale=1.d0
    e_units='au'
   elseif(units_choice==2)then
    escale=htoev
    e_units='eV'
   elseif(units_choice==3)then
    escale=htokcal
    e_units='kcal'
   elseif(units_choice==4)then
    escale=1.d0/dble(nbasis)
    e_units='au/atom'
   elseif(units_choice==5)then
    escale=htoev/dble(nbasis)
    e_units='eV/atom'
   else
    escale=htokcal/dble(nbasis)
    e_units='kcal/atom'
   endif
  else
! Default.
   do
    write(6,*)'Please select units for your energy data.'
    write(6,*)'Choose one of: (1) au; (2) eV; (3) kcal.'
    read(5,*,iostat=ierr)units_choice
    if(ierr/=0)units_choice=-1
    if(units_choice>=1.and.units_choice<=3)then
     exit
    else
     write(6,*)'Please try again.  Choose a number between 1 and 3.'
     write(6,*)
    endif
   enddo ! choice loop
   if(units_choice==1)then
    escale=1.d0
    e_units='au'
   elseif(units_choice==2)then
    escale=htoev
    e_units='eV'
   else
    escale=htokcal
    e_units='kcal'
   endif
   write(6,*)'For finite systems, energies are quoted for the whole system.'
   write(6,*)'For real crystals, energies are quoted per primitive cell.'
   write(6,*)'For electron(-hole) systems, energies are quoted per particle.'
  endif
  write(6,*)

! Rescale energy data.  (Easier just to rescale the raw data than to rescale
! each result quoted.)
  if(escale/=1.d0)then
   if(tag_energy>0)data_array(:,tag_energy)=data_array(:,tag_energy)*escale
   if(tag_etotalt>0)data_array(:,tag_etotalt)=data_array(:,tag_etotalt)*escale
   if(tag_esqr>0)data_array(:,tag_esqr)=data_array(:,tag_esqr)*escale*escale
   if(tag_popavgsqr>0)data_array(:,tag_popavgsqr)=data_array(:,tag_popavgsqr) &
    &*escale*escale
   if(tag_K>0)data_array(:,tag_K)=data_array(:,tag_K)*escale
   if(tag_T>0)data_array(:,tag_T)=data_array(:,tag_T)*escale
   if(tag_fisq>0)data_array(:,tag_fisq)=data_array(:,tag_fisq)*escale
   if(tag_Ewald>0)data_array(:,tag_Ewald)=data_array(:,tag_Ewald)*escale
   if(tag_local>0)data_array(:,tag_local)=data_array(:,tag_local)*escale
   if(tag_nonlocal>0)data_array(:,tag_nonlocal)=data_array(:,tag_nonlocal) &
    &*escale
   if(tag_short>0)data_array(:,tag_short)=data_array(:,tag_short)*escale
   if(tag_long>0)data_array(:,tag_long)=data_array(:,tag_long)*escale
   if(tag_cppei>0)data_array(:,tag_cppei)=data_array(:,tag_cppei)*escale
   if(tag_cppe>0)data_array(:,tag_cppe)=data_array(:,tag_cppe)*escale
   if(tag_cppee>0)data_array(:,tag_cppee)=data_array(:,tag_cppee)*escale
   if(tag_masspol>0)data_array(:,tag_masspol)=data_array(:,tag_masspol)*escale
   if(forces)then
    do iion=1,nitot_forces
     do iaxis=1,3
      do item=1,22
       if(tag_forces(item,iaxis,iion)/=-1)then
        data_array(:,tag_forces(item,iaxis,iion))=&
         &data_array(:,tag_forces(item,iaxis,iion))*escale
       endif
      enddo
     enddo
    enddo
   endif ! if forces
   if(tag_future0>0)data_array(:,tag_future0)=data_array(:,tag_future0)*escale
   if(tag_future1>0)data_array(:,tag_future1)=data_array(:,tag_future1)*escale
   if(tag_future2>0)data_array(:,tag_future2)=data_array(:,tag_future2)*escale
   if(tag_future3>0)data_array(:,tag_future3)=data_array(:,tag_future3)*escale
   if(tag_future4>0)data_array(:,tag_future4)=data_array(:,tag_future4)*escale
   if(tag_future5>0)data_array(:,tag_future5)=data_array(:,tag_future5)*escale
   if(tag_future6>0)data_array(:,tag_future6)=data_array(:,tag_future6)*escale
   if(tag_future7>0)data_array(:,tag_future7)=data_array(:,tag_future7)*escale
   if(tag_future8>0)data_array(:,tag_future8)=data_array(:,tag_future8)*escale
   if(tag_future9>0)data_array(:,tag_future9)=data_array(:,tag_future9)*escale
   if(tag_future10>0)data_array(:,tag_future10)=data_array(:,tag_future10)*&
    &escale
   if(tag_massvel>0)data_array(:,tag_massvel)=data_array(:,tag_massvel)*escale
   if(tag_darwinen>0)data_array(:,tag_darwinen)=data_array(:,tag_darwinen) &
    &*escale
   if(tag_darwinee>0)data_array(:,tag_darwinee)=data_array(:,tag_darwinee) &
    &*escale
   if(tag_retard>0)data_array(:,tag_retard)=data_array(:,tag_retard)*escale
   if(tag_eref>0)data_array(:,tag_eref)=data_array(:,tag_eref)*escale
   if(tag_ebest>0)data_array(:,tag_ebest)=data_array(:,tag_ebest)*escale
   constant_energy=constant_energy*escale
  endif ! Data needs rescaling.

 END SUBROUTINE check_data


 SUBROUTINE compute_stats
!--------------------------------------------------------------------------!
! In this subroutine, the various columns of data are subjected to various !
! statistical analyses.                                                    !
!--------------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER ierr,block_length,Nstudy,startline,nthird,nthirdstart, &
   &ialloc,nthirdstop,i
  DOUBLE PRECISION av,av_energy,std_err,std_err_energy,delta_std_err,var, &
   &max_val,min_val,skew,kurt,corr_tau,corr_tau_err,sqrt_tau,err_sqrt_tau, &
   &raw_var,raw_var_err,pop_var,pop_var_err
  DOUBLE PRECISION,ALLOCATABLE :: temp_data(:)

  Nstudy=Nlines-Nskip
  startline=Nskip+1

! Write out some information about the more important DMC simulation params.
! Do this first, so that important data appears at end of output.
  if(tag_nconf>0)then

   write(6,*)'ANALYSIS OF CONFIGURATION POPULATION'
   write(6,*)'===================================='
   call compute_stats_unweighted(.false.,Nstudy,data_array(startline:Nlines, &
    &tag_nconf),av,var,skew,kurt,max_val,min_val)
   write(6,*)'Minimum population : ',min_val
   write(6,*)'   Mean population : ',av
   write(6,*)'Maximum population : ',max_val
   write(6,*)'         Std error : ',sqrt(var/dble(Nstudy))
   if(av-min_val>0.25d0*av.or.max_val-av>0.25d0*av)write(6,*) &
    &'Warning: Population fluctuated by more than 25% of mean.'
   write(6,*)

  endif ! Config population data present.

  if(tag_acc>0)then

   write(6,*)'ANALYSIS OF ACCEPTANCE RATIO'
   write(6,*)'============================'
   call compute_stats_unweighted(.false.,Nstudy,data_array(startline:Nlines, &
    &tag_acc),av,var,skew,kurt,max_val,min_val)
   write(6,*)'Minimum acceptance ratio : ',min_val
   write(6,*)'   Mean acceptance ratio : ',av
   write(6,*)'Maximum acceptance ratio : ',max_val
   write(6,*)'               Std error : ',sqrt(var/dble(Nstudy))
   write(6,*)

  endif ! Acceptance-ratio data present.

  if(tag_teff>0)then

   write(6,*)'ANALYSIS OF EFFECTIVE TIME STEP'
   write(6,*)'==============================='
   call compute_stats_unweighted(.false.,Nstudy,data_array(startline:Nlines, &
    &tag_teff),av,var,skew,kurt,max_val,min_val)
   write(6,*)'Minimum time step (au) : ',min_val
   write(6,*)'   Mean time step (au) : ',av
   write(6,*)'Maximum time step (au) : ',max_val
   write(6,*)'        Std error (au) : ',sqrt(var/dble(Nstudy))
   write(6,*)

  endif ! Effective time step data present.

  if(tag_energy>0)then

! Compute lots of information about the total energy data.
   write(6,*)'ANALYSIS OF TOTAL-ENERGY DATA'
   write(6,*)'============================='
   call compute_stats_unweighted(.true.,Nstudy,data_array(startline:Nlines, &
    &tag_energy),av,var,skew,kurt,max_val,min_val)
   write(6,*)'Minimum energy (',trim(e_units),') : ',min_val
   write(6,*)'   Mean energy (',trim(e_units),') : ',av
   write(6,*)'Maximum energy (',trim(e_units),') : ',max_val
   write(6,*)'      Variance (',trim(e_units),') : ',var
   write(6,*)'     Std error (',trim(e_units),') : ',sqrt(var/dble(Nstudy))
   write(6,*)repeat(' ',len_trim(e_units))//'         Skewness : ',skew
   write(6,*)repeat(' ',len_trim(e_units))//'Normal sk. fluct. : ',&
    &sqrt(6.d0/dble(Nstudy))
   write(6,*)repeat(' ',len_trim(e_units))//'         Kurtosis : ',kurt
   write(6,*)repeat(' ',len_trim(e_units))//'Normal ku. fluct. : ',&
    &sqrt(24.d0/dble(Nstudy))
   write(6,*)'(NB, the var of the energy data is not an estimate &
    &of the actual var.)'
   write(6,*)
! Analyse total energy by thirds if there is enough data.
   if(Nstudy>=6)then
    write(6,*)'ANALYSIS OF TOTAL-ENERGY DATA BY THIRDS'
    write(6,*)'======================================='
    write(6,*)'(Energy data in units of '//trim(e_units)//'.)'
    nthird=Nstudy/3
    write(6,*)'    Data range    Av energy    Std error    Maximum      &
     &Minimum'
    do i=1,3
     nthirdstart=startline+(i-1)*nthird
     nthirdstop=nthirdstart+nthird-1
     call compute_stats_unweighted(.false.,nthird, &
      &data_array(nthirdstart:nthirdstop,tag_energy),av,var,skew,kurt, &
      &max_val,min_val)
     write(6,'(" ",a16,4(" ",es12.5))')trim(i2s(nthirdstart)) &
      &//'->'//trim(i2s(nthirdstop)),av,sqrt(var/dble(nthird)),max_val,min_val
    enddo ! i
   else
    write(6,*)'Not enough data to analyse by thirds: need at least 6 points.'
   endif ! Enough data?
   write(6,*)
   write(6,*)'CORRELATION-TIME ANALYSIS OF TOTAL-ENERGY DATA'
   write(6,*)'=============================================='
   call correlation_time(Nstudy,data_array(startline:Nlines,tag_energy), &
    &corr_tau,corr_tau_err,av,var)
   if(corr_tau/=-1.d0)then
    write(6,*)'          Correlation time (steps) : ',corr_tau
    write(6,*)' Error in correlation time (steps) : ',corr_tau_err
    write(6,*)
    if(corr_tau>0.d0)then
     sqrt_tau=sqrt(corr_tau)
     err_sqrt_tau=corr_tau_err/(2*sqrt_tau)
     write(6,*)'          Error-bar factor : ',sqrt_tau
     write(6,*)' Error in error-bar factor : ',err_sqrt_tau
     write(6,*)
     if(tag_weight>0.and.use_weights)then
      call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_energy), &
       &data_array(startline:Nlines,tag_weight),1,av,std_err,delta_std_err)
     else
      call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_energy), &
      &1,av,std_err,delta_std_err)
     endif
     write(6,*)'                 Mean energy (',trim(e_units),') : ',av
     write(6,*)' Correlation-corrected error (',trim(e_units),') : ',&
      &std_err*sqrt(corr_tau)
     write(6,*)'              Error in error (',trim(e_units),') : ',&
      &sqrt((std_err*err_sqrt_tau)**2+(sqrt_tau*delta_std_err)**2)
    else
     write(6,*)'The correlation time appears to be negative.'
    endif
   else
    write(6,*)'The correlation time could not be computed.'
   endif ! corr_tau calculated.
   write(6,*)

   write(6,*)'REBLOCKING ANALYSIS OF TOTAL-ENERGY DATA'
   write(6,*)'========================================'
! Print out reblocking analysis of energy
   write(6,*)'(Energy data in units of '//trim(e_units)//'.)'
   if(tag_weight>0.and.use_weights)then
    call reblock_analysis(Nstudy,data_array(startline:Nlines,tag_energy), &
     &data_array(startline:Nlines,tag_weight))
   else
    call reblock_analysis(Nstudy,data_array(startline:Nlines,tag_energy))
   endif ! weights

  endif ! energy data available

  do
   write(6,*)'Please choose a block length for reblocking all energy &
    &components.'
   read(5,*,iostat=ierr)block_length
   if(ierr/=0)block_length=-1
   if(block_length>=1.and.block_length<=Nstudy/2)then
    exit
   else
    write(6,*)'Please try again.  Block length should be between 1 and ' &
     &//trim(i2s(Nstudy/2))//'.'
   endif
  enddo ! get block length
  write(6,*)'Chosen block length: '//trim(i2s(block_length))//'.'
  write(6,*)

! Write out the energy components with reblocked error bars.

  write(6,*)'ENERGY COMPONENTS WITH REBLOCKED ERROR BARS'
  write(6,*)'==========================================='

5 format(" ",a30,2(" ",es22.14))
10 format(" ",a30," ",es22.14)
15 format(32x,a23,a)
  write(6,15)'  Mean ('//trim(e_units)//')                 ','  Err (' &
   &//trim(e_units)//')'

  allocate(temp_data(startline:Nlines),stat=ialloc)
  if(ialloc/=0)then
   write(6,*)'Allocation problem.'
   stop
  endif

  if(tag_weight>0.and.use_weights)then

   if(tag_energy>0)then

    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_energy), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    av_energy=av ; std_err_energy=std_err
    if(trim(interaction)=='mpc'.or.trim(interaction)=='mpc_ewald')then
     write(6,5)'Total energy (using MPC) :',av,std_err
     if(trim(interaction)=='mpc_ewald'.and.tag_ewald>0.and.tag_short>0.and.&
      &tag_long>0)then
      temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
       &-data_array(startline:Nlines,tag_short)-data_array(startline:Nlines, &
       &tag_long)+data_array(startline:Nlines,tag_ewald)
      call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
       &delta_std_err)
      write(6,5)'Total energy (using Ewald) :',av,std_err
     endif ! ewald present
    else
     if(isperiodic)then
      write(6,5)'Total energy (using Ewald) :',av,std_err
     else
      write(6,5)'Total energy :',av,std_err
     endif ! periodic
     if(trim(interaction)=='ewald_mpc'.and.tag_ewald>0.and.tag_short>0.and.&
      &tag_long>0)then
      temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
       &+data_array(startline:Nlines,tag_short)+data_array(startline:Nlines, &
       &tag_long)-data_array(startline:Nlines,tag_ewald)
      call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
       &delta_std_err)
      write(6,5)'Total energy (using MPC) :',av,std_err
     endif ! MPC present
    endif ! use_mpc_energy

    if(tag_esqr>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_esqr), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Total energy squared:',av,std_err
     raw_var=av-av_energy*av_energy
     raw_var_err=sqrt(std_err**2+(av_energy*std_err_energy)**2)
     write(6,5)'Variance of total energy:',raw_var,raw_var_err
    endif

    if(tag_popavgsqr>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_popavgsqr), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Population avg of total energy squared:',av,std_err
     pop_var=av-av_energy*av_energy
     pop_var_err=sqrt(std_err**2+(av_energy*std_err_energy)**2)
     write(6,5)'Variance of population avg of energy:',pop_var,pop_var_err
    endif

    if(tag_esqr>0.and.tag_popavgsqr>0)then
     write(6,5)'Effective population size:',raw_var/pop_var,&
      &sqrt((raw_var_err/pop_var)**2+(pop_var_err*raw_var/pop_var**2)**2)
    endif

    if(tag_masspol>0.and.tag_massvel>0.and.tag_darwinen>0.and.tag_darwinee>0 &
     &.and.tag_retard>0)then
! At present, only have relativistic data for atoms.
     temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
      &+data_array(startline:Nlines,tag_masspol) &
      &+data_array(startline:Nlines,tag_massvel) &
      &+data_array(startline:Nlines,tag_darwinen) &
      &+data_array(startline:Nlines,tag_darwinee) &
      &+data_array(startline:Nlines,tag_retard)
     call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Total energy (inc rel) :',av,std_err
    endif ! rel_present

    if(tag_K>0)then
     temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
      &-data_array(startline:Nlines,tag_K)
     call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     if(trim(interaction)=='mpc'.or.trim(interaction)=='mpc_ewald')then
      write(6,5)'Total pot energy (using MPC) :',av,std_err
     else
      if(isperiodic)then
       write(6,5)'Tot pot energy (using Ewald) :',av,std_err
      else
       write(6,5)'Total potential energy :',av,std_err
      endif ! periodic
     endif ! use_mpc_energy
    endif ! K present

   endif ! energy present.

   if(tag_K>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_K), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Kinetic energy (K) :',av,std_err
   endif ! K present
   if(tag_T>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_T), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Kinetic energy (T) :',av,std_err
   endif ! T present
   if(tag_fisq>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_fisq), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Kinetic energy (FISQ) :',av,std_err
   elseif(tag_K>0.and.tag_T>0)then
    temp_data(startline:Nlines)=2.d0*data_array(startline:Nlines,tag_T) &
     &-data_array(startline:Nlines,tag_K)
     call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Kinetic energy (FISQ) :',av,std_err
   endif ! FISQ present.
   if(tag_ewald>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_ewald), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    if(isperiodic)then
     write(6,5)'Ewald interaction :',av,std_err
    else
     write(6,5)'Coulomb interaction :',av,std_err
    endif ! periodic
   endif ! Ewald present.
   if(tag_local>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_local), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Local e-i energy :',av,std_err
   endif ! local present
   if(tag_nonlocal>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines, &
     &tag_nonlocal),data_array(startline:Nlines,tag_weight),block_length, &
     &av,std_err,delta_std_err)
    write(6,5)'Nonlocal e-i energy :',av,std_err
   endif ! nonlocal e-i pot present

  else

   if(tag_energy>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_energy), &
     &block_length,av,std_err,delta_std_err)
    av_energy=av ; std_err_energy=std_err
    if(trim(interaction)=='mpc'.or.trim(interaction)=='mpc_ewald')then
     write(6,5)'Total energy (using MPC) :',av,std_err
     if(trim(interaction)=='mpc_ewald'.and.tag_ewald>0.and.tag_long>0.and.&
      &tag_short>0)then
      temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
       &-data_array(startline:Nlines,tag_short)-data_array(startline:Nlines, &
       &tag_long)+data_array(startline:Nlines,tag_ewald)
      call reblock_unweighted(Nstudy,temp_data(startline:Nlines), &
       &block_length,av,std_err,delta_std_err)
      write(6,5)'Total energy (using Ewald) :',av,std_err
     endif ! Ewald present.
    else
     if(isperiodic)then
      write(6,5)'Total energy (using Ewald) :',av,std_err
     else
      write(6,5)'Total energy :',av,std_err
     endif ! periodic
     if(trim(interaction)=='ewald_mpc'.and.tag_ewald>0.and.tag_long>0.and.&
      &tag_short>0)then
      temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
       &+data_array(startline:Nlines,tag_short)+data_array(startline:Nlines, &
       &tag_long)-data_array(startline:Nlines,tag_ewald)
      call reblock_unweighted(Nstudy,temp_data(startline:Nlines), &
       &block_length,av,std_err,delta_std_err)
      write(6,5)'Total energy (using MPC) :',av,std_err
     endif ! MPC present.
    endif ! use_mpc_energy
    if(tag_esqr>0)then
     call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_esqr), &
      &block_length,av,std_err,delta_std_err)
     write(6,5)'Total energy squared:',av,std_err
     write(6,10)'Raw variance of total energy:',&
      &(av-av_energy*av_energy)*dble(nparts_per_simcell)
    endif
    if(tag_masspol>0.and.tag_massvel>0.and.tag_darwinen>0.and.tag_darwinee>0 &
     &.and.tag_retard>0)then
! At present, only have relativistic data for atoms.
     temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
      &+data_array(startline:Nlines,tag_masspol) &
      &+data_array(startline:Nlines,tag_massvel) &
      &+data_array(startline:Nlines,tag_darwinen) &
      &+data_array(startline:Nlines,tag_darwinee) &
      &+data_array(startline:Nlines,tag_retard)
     call reblock_unweighted(Nstudy,temp_data(startline:Nlines),block_length, &
      &av,std_err,delta_std_err)
     write(6,5)'Total energy (inc rel) :',av,std_err
    endif ! rel_present
    if(tag_K>0)then
     temp_data(startline:Nlines)=data_array(startline:Nlines,tag_energy) &
      &-data_array(startline:Nlines,tag_K)
     call reblock_unweighted(Nstudy,temp_data(startline:Nlines),block_length, &
      &av,std_err,delta_std_err)
     if(trim(interaction)=='mpc'.or.trim(interaction)=='mpc_ewald')then
      write(6,5)'Total pot energy (using MPC) :',av,std_err
     else
      if(isperiodic)then
       write(6,5)'Tot pot energy (using Ewald) :',av,std_err
      else
       write(6,5)'Total potential energy :',av,std_err
      endif ! periodic
     endif ! use_mpc_energy
    endif ! K present
   endif ! Energy data present.
   if(tag_K>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_K), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Kinetic energy (K) :',av,std_err
   endif ! K present.
   if(tag_T>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_T), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Kinetic energy (T) :',av,std_err
   endif ! T present
   if(tag_fisq>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_fisq), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Kinetic energy (FISQ) :',av,std_err
   elseif(tag_K>0.and.tag_T>0)then
    temp_data(startline:Nlines)=2.d0*data_array(startline:Nlines,tag_T) &
     &-data_array(startline:Nlines,tag_K)
    call reblock_unweighted(Nstudy,temp_data(startline:Nlines),block_length, &
     &av,std_err,delta_std_err)
    write(6,5)'Kinetic energy (FISQ) :',av,std_err
   endif ! K & T present.
   if(tag_ewald>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_ewald), &
     &block_length,av,std_err,delta_std_err)
    if(isperiodic)then
     write(6,5)'Ewald interaction :',av,std_err
    else
     write(6,5)'Coulomb interaction :',av,std_err
    endif ! periodic
   endif ! Ewald present.
   if(tag_local>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_local), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Local e-i energy :',av,std_err
   endif ! local present
   if(tag_nonlocal>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines, &
     &tag_nonlocal),block_length,av,std_err,delta_std_err)
    write(6,5)'Nonlocal e-i energy :',av,std_err
   endif ! nonlocal e-i pot present

  endif ! weighted

  if(tag_weight>0.and.use_weights)then
   if(tag_short>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_short), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Short-range MPC energy :',av,std_err
   endif ! short present
   if(tag_long>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_long), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Long-range MPC energy :',av,std_err
   endif ! long present.
   if(tag_short>0.and.tag_long>0)then
    temp_data(startline:Nlines)=data_array(startline:Nlines,tag_short)+ &
     &data_array(startline:Nlines,tag_long)
    call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Total MPC energy :',av,std_err
   endif ! MPC data present
   if(tag_cppei>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_cppei), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'CPP energy (e-i) :',av,std_err
   endif ! CPPEI present
   if(tag_cppe>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_cppe), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'CPP energy (e) :',av,std_err
   endif ! CPPE present
   if(tag_cppee>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_cppee), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'CPP energy (e-e) :',av,std_err
   endif ! CPP data present
   if(tag_masspol>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_masspol), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Mass-polarization energy :',av,std_err
   endif ! masspol present
   if(tag_massvel>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_massvel), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Mass-velocity energy :',av,std_err
   endif ! massvel present.
   if(tag_darwinen>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_darwinen), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Electron-nucleus Darwin :',av,std_err
   endif ! darawin e-n present.
   if(tag_darwinee>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_darwinee), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Electron-electron Darwin :',av,std_err
   endif ! darwin e-e present.
   if(tag_retard>0)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_retard), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Retardation term :',av,std_err
   endif ! Duh.
   if(tag_masspol>0.and.tag_massvel>0.and.tag_darwinen>0.and.tag_darwinee>0 &
    &.and.tag_retard>0)then
    temp_data(startline:Nlines)=data_array(startline:Nlines,tag_masspol) &
     &+data_array(startline:Nlines,tag_massvel) &
     &+data_array(startline:Nlines,tag_darwinen) &
     &+data_array(startline:Nlines,tag_darwinee) &
     &+data_array(startline:Nlines,tag_retard)
    call reblock_weighted(Nstudy,temp_data(startline:Nlines), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Total rel correction :',av,std_err
   endif ! Relativistic data present

  else

   if(tag_short>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_short), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Short-range MPC energy :',av,std_err
   endif ! short present.
   if(tag_long>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_long), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Long-range MPC energy :',av,std_err
   endif ! long present.
   if(tag_short>0.and.tag_long>0)then
    temp_data(startline:Nlines)=data_array(startline:Nlines,tag_short) &
     &+data_array(startline:Nlines,tag_long)
    call reblock_unweighted(Nstudy,temp_data(startline:Nlines),block_length, &
     &av,std_err,delta_std_err)
    write(6,5)'Total MPC energy :',av,std_err
   endif ! MPC data present
   if(tag_cppei>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_cppei), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'CPP energy (e-i) :',av,std_err
   endif ! CPPEI present
   if(tag_cppe>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_cppe), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'CPP energy (e) :',av,std_err
   endif ! CPPE present.
   if(tag_cppee>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_cppee), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'CPP energy (e-e) :',av,std_err
   endif ! CPPEE data present
   if(tag_masspol>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_masspol), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Mass-polarization energy :',av,std_err
   endif ! masspol present.
   if(tag_massvel>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_massvel), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Mass-velocity energy :',av,std_err
   endif ! massvel present.
   if(tag_darwinen>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_darwinen), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Electron-nucleus Darwin :',av,std_err
   endif ! darwinen present.
   if(tag_darwinee>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_darwinee), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Electron-electron Darwin :',av,std_err
   endif ! darwinee present.
   if(tag_retard>0)then
    call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_retard), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Retardation term :',av,std_err
   endif ! Duh.
   if(tag_masspol>0.and.tag_massvel>0.and.tag_darwinen>0.and.tag_darwinee>0 &
    &.and.tag_retard>0)then
    temp_data(startline:Nlines)=data_array(startline:Nlines,tag_masspol) &
     &+data_array(startline:Nlines,tag_massvel) &
     &+data_array(startline:Nlines,tag_darwinen) &
     &+data_array(startline:Nlines,tag_darwinee) &
     &+data_array(startline:Nlines,tag_retard)
    call reblock_unweighted(Nstudy,temp_data(startline:Nlines), &
     &block_length,av,std_err,delta_std_err)
    write(6,5)'Total rel correction :',av,std_err
   endif ! Relativistic data present

  endif ! weighted.

  deallocate(temp_data)

  if(constant_energy/=0.d0)write(6,10)'Constant energy :',constant_energy
  write(6,*)

! Calculate and write out forces
  if(forces)call construct_write_forces(startline,Nstudy,block_length)

! Write out future-walking estimates
  if(tag_future1>0)then
   write(6,*)'FUTURE-WALKING ESTIMATES WITH REBLOCKED ERROR BAR'
   write(6,*)'================================================='
   write(6,*)'Future-walking estimates of the observable pureitems(1) &
    &from the dmc_main'
   write(6,*)'routine. Temporarily, pureitems(1) is the Hellmann-Feynman &
    &forces in the'
   write(6,*)'x-direction of the 1st atom as ordered in the gwfn.data &
    &file. To estimate'
   write(6,*)'a different observable, alter the assignment after line:'
   write(6,*)"'Change next line when future-walking estimates are required'"
   write(6,*)'                                 Mean (au)              Err (au)'

   if(tag_weight>0.and.use_weights)then
    if(tag_future0>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future0), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (1st) :',av,std_err
    endif
    if(tag_future1>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future1), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (2nd) :',av,std_err
    endif
    if(tag_future2>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future2), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (3rd) :',av,std_err
    endif
    if(tag_future3>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future3), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (4th) :',av,std_err
    endif
    if(tag_future4>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future4), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (5th) :',av,std_err
    endif
    if(tag_future5>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future5), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (6th) :',av,std_err
    endif
    if(tag_future6>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future6), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (7th) :',av,std_err
    endif
     if(tag_future7>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future7), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator  (8th) :',av,std_err
    endif
    if(tag_future8>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future8), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
       &delta_std_err)
     write(6,5)'FW Estimator  (9th) :',av,std_err
    endif
    if(tag_future9>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future9), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator (10th) :',av,std_err
    endif
    if(tag_future10>0) then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_future10), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'FW Estimator (11th) :',av,std_err
    endif
    write(6,*)
    write(6,*)'The 1st estimator corresponds to future-walking projection &
     &time T=0 1/Ha,'
    write(6,*)'the 2nd to T=0.5 1/Ha, the 3rd to T=1 1/Ha..., and the &
     &11th to T=10 1/Ha.'
    write(6,*)

   endif ! tag_weight

  endif ! future

  if(tag_dipole1>0.or.tag_dipole2>0.or.tag_dipole3>0.or.tag_dipole_sq>0)then
   write(6,*)'ELECTRIC DIPOLE MOMENT WITH REBLOCKED ERROR BARS'
   write(6,*)'================================================'
   write(6,15)'  Mean (au)                 ','  Err (au)'

   if(tag_weight>0.and.use_weights)then
    if(tag_dipole1>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_dipole1), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Dipole moment (x cpt) :',av,std_err
    endif ! tag_dipole1
    if(tag_dipole2>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_dipole2), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Dipole moment (y cpt) :',av,std_err
    endif ! tag_dipole2
    if(tag_dipole3>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_dipole3), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Dipole moment (z cpt) :',av,std_err
    endif ! tag_dipole3
    if(tag_dipole_sq>0)then
     call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_dipole_sq), &
      &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
      &delta_std_err)
     write(6,5)'Dipole moment squared :',av,std_err
    endif ! tag_dipole_sq
   else
    if(tag_dipole1>0)then
     call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_dipole1), &
      &block_length,av,std_err,delta_std_err)
     write(6,5)'Dipole moment (x cpt) :',av,std_err
    endif ! tag_dipole1
    if(tag_dipole2>0)then
     call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_dipole2), &
      &block_length,av,std_err,delta_std_err)
     write(6,5)'Dipole moment (y cpt) :',av,std_err
    endif ! tag_dipole2
    if(tag_dipole3>0)then
     call reblock_unweighted(Nstudy,data_array(startline:Nlines,tag_dipole3), &
      &block_length,av,std_err,delta_std_err)
     write(6,5)'Dipole moment (z cpt) :',av,std_err
    endif ! tag_dipole3
    if(tag_dipole_sq>0)then
     call reblock_unweighted(Nstudy,data_array(startline:Nlines, &
      &tag_dipole_sq),block_length,av,std_err,delta_std_err)
     write(6,5)'Dipole moment squared :',av,std_err
    endif ! tag_dipole_sq
   endif ! weighted.
   write(6,*)
  endif ! Dipole moment

  if(tag_contact_den>0)then
   write(6,*)'CONTACT DENSITY'
   write(6,*)'==============='
   write(6,15)'  Mean (au)                 ','  Err (au)'
   if(tag_weight>0.and.use_weights)then
    call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_contact_den), &
     &data_array(startline:Nlines,tag_weight),block_length,av,std_err, &
     &delta_std_err)
    write(6,5)'Elec-pos contact density :',av,std_err
   else
    call reblock_unweighted(Nstudy,data_array(startline:Nlines, &
     &tag_contact_den),block_length,av,std_err,delta_std_err)
    write(6,5)'Elec-pos contact density :',av,std_err
   endif ! weighted.
   write(6,*)
  endif ! Contact density

 END SUBROUTINE compute_stats


 SUBROUTINE reblock_analysis(no_pts,data_array,weight_array)
!--------------------------------------------------------------!
! Compute the weighted average of the data, and calculate the  !
! error bar as a function of reblocking transformation number. !
!--------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: no_pts
  DOUBLE PRECISION,INTENT(in) :: data_array(no_pts)
  DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight_array(no_pts)
  INTEGER no_rtns,rtn,block_length,ierr
  DOUBLE PRECISION av,std_err,delta_std_err

  open(unit=10,file='reblock.plot',status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening reblock.plot.'
   stop
  endif

! Number of reblocking transformations
  no_rtns=floor(log(dble(no_pts))/log(2.d0))

! Initial block length
  block_length=1

! Write out results of reblocking analysis
  write(6,*)' RTN   Blk leng   Std error in mean       Error in std error'
  do rtn=0,no_rtns-1
   if(present(weight_array))then
    call reblock_weighted(no_pts,data_array,weight_array,block_length, &
     &av,std_err,delta_std_err)
   else
    call reblock_unweighted(no_pts,data_array,block_length,av,std_err, &
     &delta_std_err)
   endif ! weights present
   write(6,'(" ",i4," ",i10," ",es23.15," ",es23.15)')rtn,block_length, &
    &std_err,delta_std_err
   write(10,*)rtn,std_err,delta_std_err
   block_length=2*block_length
  enddo ! rtn
  write(6,*)

  write(6,*)'Reblocked error bar against reblocking transformation number &
   &(RTN) has been'
  write(6,*)'written to reblock.plot.  Please use "plot_reblock" to view &
   &these data.'
  write(6,*)
  close(10)

 END SUBROUTINE reblock_analysis


 SUBROUTINE reblock_forces_analysis(no_pts,data_array,plotname,weight_array)
!--------------------------------------------------------------!
! Compute the weighted average of the data, and calculate the  !
! error bar as a function of reblocking transformation number. !
! This routine is an extension of routine reblock_analysis     !
! to allow specifying the name of the file to be written out.  !
!                                                              !
! AB   11.2007                                                 !
!--------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: no_pts
  DOUBLE PRECISION,INTENT(in) :: data_array(no_pts)
  DOUBLE PRECISION,INTENT(in),OPTIONAL :: weight_array(no_pts)
  CHARACTER(20),INTENT(in) :: plotname
  INTEGER no_rtns,rtn,block_length,ierr
  DOUBLE PRECISION av,std_err,delta_std_err

  open(unit=10,file=plotname ,status='replace',iostat=ierr)
  if(ierr/=0)then
   write(6,*)'Problem opening plotname.'
   stop
  endif

! Number of reblocking transformations
  no_rtns=floor(log(dble(no_pts))/log(2.d0))

! Initial block length
  block_length=1

! Write out results of reblocking analysis
  write(6,*)' RTN   Blk leng   Std error in mean       Error in std error'
  do rtn=0,no_rtns-1
   if(present(weight_array))then
    call reblock_weighted(no_pts,data_array,weight_array,block_length, &
     &av,std_err,delta_std_err)
   else
    call reblock_unweighted(no_pts,data_array,block_length,av,std_err, &
     &delta_std_err)
   endif ! weights present
   write(6,'(" ",i4," ",i10," ",es23.15," ",es23.15)')rtn,block_length, &
    &std_err,delta_std_err
   write(10,*)rtn,std_err,delta_std_err
   block_length=2*block_length
  enddo ! rtn
  write(6,*)

  close(10)

 END SUBROUTINE reblock_forces_analysis


 CHARACTER(12) FUNCTION i2s(n)
!-----------------------------------------------------------------------!
! I2S                                                                   !
! ===                                                                   !
! Convert integers to left justified strings that can be printed in the !
! middle of a sentence without introducing large amounts of white space.!
!                                                                       !
! Calling routine is intended to include something like:                !
! USE utilities                                                         !
! INTEGER i                                                             !
! i=12                                                                  !
! write(6,*)'Integer number ',trim(i2s(i)),' with words at the end.'    !
!-----------------------------------------------------------------------!
  IMPLICIT NONE
  INTEGER,INTENT(in) :: n
  INTEGER i,j
  INTEGER,PARAMETER :: ichar0=ichar('0')
  i2s=''
  i=abs(n)
  do j=len(i2s),1,-1
   i2s(j:j)=achar(ichar0+mod(i,10))
   i=i/10 ; if(i==0)exit
  enddo ! j
  if(n<0)then
   i2s='-'//adjustl(i2s)
  else
   i2s=adjustl(i2s)
  endif ! n<0
 END FUNCTION i2s


 SUBROUTINE generate_tag_forces(datastring,i)
!------------------------------------------------------------------!
! This subroutine generates tags for forces from tags read in from !
! the .hist file. It also determines the number of atoms and axes  !
! for which forces data are available.                             !
!                                                                  !
! AB 11.2007                                                       !
!------------------------------------------------------------------!
  IMPLICIT NONE
  CHARACTER(72),INTENT(in) :: datastring
  INTEGER,INTENT(in) :: i
  CHARACTER(1) :: axis(3)=(/'X','Y','Z'/)
  CHARACTER(1) :: atem(22)=(/'A','B','C','D','E','F','G','H','I','J','K',&
   &'L','M','N','O','P','Q','R','S','T','U','V'/)
  CHARACTER(1) iaxis_char,item_char
  INTEGER iaxis,item,iion,ion_tmp

  read(datastring(6:),*)ion_tmp
  do iion=1,nitot_max_forces
   if(ion_tmp==iion)then
! Generate the number of atoms for which forces are calculated.
    if(iion>nitot_forces)nitot_forces=iion
    do iaxis=1,3
     iaxis_char=axis(iaxis)
     if(datastring(5:5)==iaxis_char)then
! Generate the number of axis for which forces are calculated.
      if(iaxis>naxis_forces)naxis_forces=iaxis
      do item=1,22
       item_char=atem(item)
       if(datastring(4:4)==item_char)then
! Generate the number of items.
        if(item>nitem_forces)nitem_forces=item
        call check_tag_free(tag_forces(item,iaxis,iion))
! Generate label for forces
        tag_forces(item,iaxis,iion)=i
       endif ! item_char
      enddo ! item
     endif ! datastring=iaxis_char
    enddo ! iaxis
   endif ! iion
  enddo ! iion

 END SUBROUTINE generate_tag_forces


 SUBROUTINE construct_write_forces(startline,Nstudy,block_length)
!----------------------------------------------------------------!
! This routine calculates VMC/DMC forces from available data and !
! performs a reblocking analysis.                                !
!                                                                !
! AB 11.2007                                                     !
!----------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: startline,Nstudy,block_length
 CHARACTER(20) plotname
 INTEGER i,n,ialloc,nthird,nthirdstart,nthirdstop,ierr
 DOUBLE PRECISION etot,std_err,delta_std_err,etot_dmc,etot_dmc_SE
 DOUBLE PRECISION av,var,skew,kurt,max_val,min_val
 LOGICAL forces_reblock,ltemp

! Do we want to reblock forces?
   forces_reblock=.false.
   do
    write(6,*)'Forces data are detected. When you like to reblock the forces da&
     &ta with the'
    write(6,*)'same block length as the total energy, choose F. When you like t&
     &o investigate'
    write(6,*)'reblocked forces error bars, choose T and use gnuplot to look at&
     & error bars.'
    write(6,*)'Choose F or T:'
    read(5,*,iostat=ierr)forces_reblock
    if(ierr/=0)forces_reblock=.false.
    if((forces_reblock).or.(.not.forces_reblock))then
     exit
    else
     write(6,*)'Please try again.  Choose T or F.'
    endif
   enddo ! choice loop


  write(6,*)'FORCES COMPONENTS WITH REBLOCKED ERROR BARS'
  write(6,*)'==========================================='

!------------------ reblock VMC forces -------------------------

  if(trim(qmc_method)=='VMC')then

   allocate(forces_array(Nlines,11),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)'Allocation problem (1).'
    stop
   endif

! Need energy estimate
   call reblock_unweighted(Nstudy,data_array(startline:Nlines,&
    &tag_energy),1,etot,std_err,delta_std_err)

   inquire(file='DMC_energy',exist=ltemp)
   if(ltemp)then
    open(11,file='DMC_energy')
    read(11,*)etot_dmc,etot_dmc_SE
    close(11)
   else
    etot_dmc=0.d0
   endif

   do iion=1,nitot_forces
    do iaxis=1,naxis_forces
     write(6,*)'Forces on atom ',trim(i2s(iion)),' along axis ',&
      &trim(i2s(iaxis)),'        Mean (au)              Err (au)'

! Construct various VMC estimators for the forces
     do n=1,Nlines
! 1. Total forces (d-loc)
      forces_array(n,1)=data_array(n,tag_forces(9,iaxis,iion))&
       &-2.d0*data_array(n,tag_forces(2,iaxis,iion))&
       &+2.d0*etot*data_array(n,tag_forces(1,iaxis,iion))
!    HFT forces (d-loc)
      forces_array(n,2)=data_array(n,tag_forces(9,iaxis,iion))
      if(data_array(n,tag_forces(10,iaxis,iion))/=0.d0)then
!    HFT forces (p-loc)
       forces_array(n,3)=data_array(n,tag_forces(10,iaxis,iion))
!    HFT forces (s-loc)
       forces_array(n,4)=data_array(n,tag_forces(11,iaxis,iion))
      endif
!    Wavefunction Pulay term
      forces_array(n,5)=&
       &-2.d0*data_array(n,tag_forces(2,iaxis,iion))&
       &+2.d0*etot*data_array(n,tag_forces(1,iaxis,iion))
!    Pseudopotential Pulay term
      forces_array(n,6)=-data_array(n,tag_forces(7,iaxis,iion))&
       &+data_array(n,tag_forces(4,iaxis,iion))
! 2. Total forces zero-variance corrected (class 1,d-loc)
      forces_array(n,7)=data_array(n,tag_forces(9,iaxis,iion))& ! HFT
       &-2.d0*data_array(n,tag_forces(2,iaxis,iion))&           ! Pulay
       &+2.d0*etot*data_array(n,tag_forces(1,iaxis,iion))&      !   "
       &-data_array(n,tag_forces(6,iaxis,iion))&                ! -H Psi'
       &-data_array(n,tag_forces(7,iaxis,iion))&                !   "
       &+data_array(n,tag_forces(3,iaxis,iion))&                !   "
       &+data_array(n,tag_forces(4,iaxis,iion))                 !   "
!    Zero-variance term
      forces_array(n,8)=&
       &-data_array(n,tag_forces(6,iaxis,iion))&                !- H Psi'
       &-data_array(n,tag_forces(7,iaxis,iion))&                !   "
       &+data_array(n,tag_forces(3,iaxis,iion))&                !   "
       &+data_array(n,tag_forces(4,iaxis,iion))                 !   "
                                                      ! E_l Psi' term cancelled
! VMC nodal term, added to Total Forces (purHFT,purNT,d-loc)
      if(etot_dmc/=0.d0)then
       forces_array(n,9)=&                                       !
        &-data_array(n,tag_forces(6,iaxis,iion))&                ! -H Psi'
        &-data_array(n,tag_forces(7,iaxis,iion))&                !   "
        &-data_array(n,tag_forces(2,iaxis,iion))&                !   "
        &+data_array(n,tag_forces(3,iaxis,iion))&                !   "
        &+data_array(n,tag_forces(4,iaxis,iion))&                !   "
        &+etot_dmc*data_array(n,tag_forces(1,iaxis,iion))        !  E*Psi'
      else
       forces_array(n,9)=0.d0
      endif
     enddo ! Nlines

! Write out the various estimators for the forces
! 1. Total forces
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,1),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'Total Force(dloc) :',av,std_err
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,2),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'HFT Force(dloc) :',av,std_err
     if(data_array(n,tag_forces(10,iaxis,iion))/=0.d0)then
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,3),&
       &block_length,av,std_err,delta_std_err)
      write(6,9) 'HFT Force(ploc) :',av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,4),&
       &block_length,av,std_err,delta_std_err)
      write(6,9) 'HFT Force(sloc) :',av,std_err
     endif
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,5),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'Wavefunction Pulay term :',av,std_err
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,6),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'Pseudopotential Pulay term :',av,std_err

! 2. Total forces zero-variance corrected
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,7),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'Total Force+ZV(dloc) :',av,std_err
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,8),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'Zero-variance term :',av,std_err
     call reblock_unweighted(Nstudy,forces_array(startline:Nlines,9),&
      &block_length,av,std_err,delta_std_err)
     write(6,9) 'VMC NT(add to last DMC est) :',av,std_err

     if(forces_reblock)then
      plotname='forces'//'.tot.'//trim(i2s(iaxis))//trim(i2s(iion))//&
       &'.plot'
      call reblock_forces_analysis(Nstudy,forces_array(startline:&
       &Nlines,1),plotname)
      plotname='forces'//'.HFT.'//trim(i2s(iaxis))//trim(i2s(iion))//&
       &'.plot'
      call reblock_forces_analysis(Nstudy,forces_array(startline:&
       &Nlines,2),plotname)
      plotname='forces'//'.totZV.'//trim(i2s(iaxis))//trim(i2s(iion))//&
       &'.plot'
      call reblock_forces_analysis(Nstudy,forces_array(startline:&
       &Nlines,7),plotname)
      plotname='forces'//'.vmcNT.'//trim(i2s(iaxis))//trim(i2s(iion))//&
       &'.plot'
      call reblock_forces_analysis(Nstudy,forces_array(startline:&
       &Nlines,9),plotname)
! Analyse total forces by thirds if there is enough data.
      if(Nstudy>=6)then
       write(6,*)
       write(6,*)'ANALYSIS OF TOTAL FORCES DATA BY THIRDS'
       nthird=Nstudy/3
       write(6,*)'    Data range    Av forces    Variance     Maximum      &
        &Minimum'
       do i=1,3
        nthirdstart=startline+(i-1)*nthird
        nthirdstop=nthirdstart+nthird-1
        call compute_stats_unweighted(.false.,nthird, &
         &forces_array(nthirdstart:nthirdstop,1),av,var,skew,kurt, &
         &max_val,min_val)
        write(6,'(" ",a16,4(" ",es12.5))')trim(i2s(nthirdstart)) &
         &//'->'//trim(i2s(nthirdstop)),av,var,max_val,min_val
       enddo ! i
      else
       write(6,*)'Not enough data to analyse by thirds: need at least 6 points.'
      endif ! Enough data?
     endif ! reblock
     write(6,*)

    enddo ! iaxis
   enddo ! iion
   write(6,*) "The last estimator 'VMC NT' is zero (not used), unless a 'DMC_en&
    &ergy' file is "
   write(6,*) "provided in the working directory during the reblocking process &
    &which contains"
   write(6,*) "the DMC energy followed by its error bar. After this estimator V&
    &MC NT is "
   write(6,*) "calculated by the reblocking routine, please add it by hand to t&
    &he 'Total"
   write(6,*) "Force(purHFT/NT,dloc)' estimator to obtain another force estimat&
    &or."
   write(6,*)
  endif ! VMC

!------------------- reblock DMC forces ---------------------

  if(trim(qmc_method)=='DMC')then
   allocate(forces_array(Nlines,10),stat=ialloc)
   if(ialloc/=0)then
    write(6,*)'Allocation problem (1).'
    stop
   endif

! Need energy estimate
   call reblock_weighted(Nstudy,data_array(startline:Nlines,tag_energy),&
    &data_array(startline:Nlines,tag_weight),1,etot,std_err,delta_std_err)

! Construct various DMC estimators for the forces
   do iion=1,nitot_forces
    do iaxis=1,naxis_forces
     write(6,*)'Forces on atom ',trim(i2s(iion)),' along axis ',&
      &trim(i2s(iaxis)),'        Mean (au)              Err (au)'

     do n=1,Nlines
! 1. Pure total forces (mixNT,dloc)
      forces_array(n,1)=data_array(n,tag_forces(20,iaxis,iion))& ! HFT
       &-data_array(n,tag_forces(18,iaxis,iion))&                ! PPT
       &+data_array(n,tag_forces(15,iaxis,iion))&                !  "
       &-2.d0*(data_array(n,tag_forces(6,iaxis,iion))&           ! NT
       &      +data_array(n,tag_forces(8,iaxis,iion))&           ! "
       &      +data_array(n,tag_forces(2,iaxis,iion))&           ! "
       &      -data_array(n,tag_forces(3,iaxis,iion))&           ! "
       &      -data_array(n,tag_forces(5,iaxis,iion))&           ! "
       & -etot*data_array(n,tag_forces(1,iaxis,iion)))           ! "
!    Pure total forces (purNT,dloc)
      forces_array(n,2)=data_array(n,tag_forces(20,iaxis,iion))& ! HFT
       &-data_array(n,tag_forces(18,iaxis,iion))&                ! PPT
       &+data_array(n,tag_forces(15,iaxis,iion))&                !  "
       &-(data_array(n,tag_forces(17,iaxis,iion))&           ! NT
       & +data_array(n,tag_forces(19,iaxis,iion))&           ! "
       & +data_array(n,tag_forces(13,iaxis,iion))&           ! "
       & -data_array(n,tag_forces(14,iaxis,iion))&           ! "
       & -data_array(n,tag_forces(16,iaxis,iion))&           ! "
       & -etot*data_array(n,tag_forces(12,iaxis,iion)))           ! "
!    Pure HFT forces (d-loc)
      forces_array(n,3)=data_array(n,tag_forces(20,iaxis,iion))
      if(data_array(n,tag_forces(10,iaxis,iion))/=0.d0)then
!    Pure HFT forces (p-loc)
       forces_array(n,4)=data_array(n,tag_forces(21,iaxis,iion))
!    Pure HFT forces (s-loc)
       forces_array(n,5)=data_array(n,tag_forces(22,iaxis,iion))
      endif
!    Nodal term (mix)
      forces_array(n,6)=&
       &-2.d0*(data_array(n,tag_forces(6,iaxis,iion))&           ! NT
       &      +data_array(n,tag_forces(8,iaxis,iion))&           ! "
       &      +data_array(n,tag_forces(2,iaxis,iion))&           ! "
       &      -data_array(n,tag_forces(3,iaxis,iion))&           ! "
       &      -data_array(n,tag_forces(5,iaxis,iion))&           ! "
       & -etot*data_array(n,tag_forces(1,iaxis,iion)))           ! "
!    Nodal term (pur)
      forces_array(n,7)=&
       &-(data_array(n,tag_forces(17,iaxis,iion))&          ! NT
       & +data_array(n,tag_forces(19,iaxis,iion))&          ! "
       & +data_array(n,tag_forces(13,iaxis,iion))&          ! "
       & -data_array(n,tag_forces(14,iaxis,iion))&          ! "
       & -data_array(n,tag_forces(16,iaxis,iion))&          ! "
       & -etot*data_array(n,tag_forces(12,iaxis,iion)))          ! "
!    Pseudopotential Pulay term: Psi^(-1)WPsi'-Psi^(-1)WPsi Psi'/Psi
      forces_array(n,8)=-data_array(n,tag_forces(18,iaxis,iion))&
       &+data_array(n,tag_forces(15,iaxis,iion))
! 2. Mixed total forces (d-loc)
      forces_array(n,9)=data_array(n,tag_forces(9,iaxis,iion))& ! HFT
       &      -data_array(n,tag_forces(7,iaxis,iion))&          ! PPT
       &      +data_array(n,tag_forces(4,iaxis,iion))&          !  "
       &      -data_array(n,tag_forces(6,iaxis,iion))&          ! NT
       &      -data_array(n,tag_forces(8,iaxis,iion))&          ! "
       & -2.d0*data_array(n,tag_forces(2,iaxis,iion))&          ! "
       &      +data_array(n,tag_forces(3,iaxis,iion))&          ! "
       &      +data_array(n,tag_forces(5,iaxis,iion))&          ! "
       & +2.d0*etot*data_array(n,tag_forces(1,iaxis,iion))      ! "

!    Mixed HFT forces (d-loc)
      forces_array(n,10)=data_array(n,tag_forces(9,iaxis,iion))
     enddo ! Nlines

! Write out forces
     if(tag_weight>0.and.use_weights)then
! 1. Pure total forces
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,1),&
       &data_array(startline:Nlines,tag_weight),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Total Force(purHFT,mixNT,dloc) :',av,std_err
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,2),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'Total Force(purHFT,purNT,dloc) :',av,std_err
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,3),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'HFT Force(pur,dloc) :',av,std_err
      if(data_array(n,tag_forces(10,iaxis,iion))/=0.d0)then
       call reblock_weighted(Nstudy,forces_array(startline:Nlines,4),&
        &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
        &delta_std_err)
       write(6,9)'HFT Force(pur,ploc) :',av,std_err
       call reblock_weighted(Nstudy,forces_array(startline:Nlines,5),&
        &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
        &delta_std_err)
       write(6,9)'HFT Force(pur,sloc) :',av,std_err
      endif
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,6),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'Nodal Term(mix) :', av,std_err
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,7),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'Nodal Term(pur) :', av,std_err
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,8),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'Pseudopot. Pulay Term(pur) :',av,std_err
! Mixed total forces
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,9),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'Total Force(mix,dloc) :',av,std_err
      call reblock_weighted(Nstudy,forces_array(startline:Nlines,10),&
       &data_array(startline:Nlines,tag_weight),block_length,av,std_err,&
       &delta_std_err)
      write(6,9)'HFT Force(mix,dloc) :',av,std_err
      write(6,*)

! Reblock forces and write out into files
      if(forces_reblock)then
       plotname='forces'//'.totpur.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,1),&
        &plotname,data_array(startline:Nlines,tag_weight))
       plotname='forces'//'.HFTpur.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,3),&
        &plotname,data_array(startline:Nlines,tag_weight))
       plotname='forces'//'.totmix'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,9),&
        &plotname,data_array(startline:Nlines,tag_weight))
       plotname='forces'//'.HFT'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,10),&
        &plotname,data_array(startline:Nlines,tag_weight))
      endif ! forces_reblock

     else ! have weights
! 1. Pure total forces
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,1),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Total Force(purHFT,mixHFT,dloc) :',av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,2),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Total Force(purHFT,purHFT,dloc) :',av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,3),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'HFT Force(pur,dloc) :',av,std_err
      if(data_array(n,tag_forces(10,iaxis,iion))/=0.d0)then
       call reblock_unweighted(Nstudy,forces_array(startline:Nlines,4),&
        &block_length,av,std_err,delta_std_err)
       write(6,9)'HFT Force(pur,ploc) :',av,std_err
       call reblock_unweighted(Nstudy,forces_array(startline:Nlines,5),&
        &block_length,av,std_err,delta_std_err)
       write(6,9)'HFT Force(pur,sloc) :',av,std_err
      endif
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,6),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Nodal Term(mix) :', av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,7),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Nodal Term(pur) :',av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,8),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Pseudopot. Pulay Term(pur) :',av,std_err
! Mixed total forces
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,9),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'Total Force(mix,dloc) :',av,std_err
      call reblock_unweighted(Nstudy,forces_array(startline:Nlines,10),&
       &block_length,av,std_err,delta_std_err)
      write(6,9)'HFT Force(mix,dloc) :',av,std_err
      write(6,*)

! Reblock forces and write out into files
      if(forces_reblock)then
       plotname='forces'//'.totpur.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,1),&
        &plotname)
       plotname='forces'//'.HFTpur.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,3),&
        &plotname)
       plotname='forces'//'.totmix.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,9),&
        &plotname)
       plotname='forces'//'.HFTmix.'//trim(i2s(iaxis))//trim(i2s(iion))//&
        &'.plot'
       call reblock_forces_analysis(Nstudy,forces_array(startline:Nlines,10),&
        &plotname)
      endif
     endif ! have weights

    enddo ! naxis
   enddo ! nitot_forces

  endif ! DMC

9 format(" ",a32,2(" ",f20.14))

 END SUBROUTINE construct_write_forces


END MODULE analysis


PROGRAM analyse_qmc
!---------------------------!
! Main program starts here. !
!---------------------------!
 USE analysis, ONLY : filename,read_data,check_data,compute_stats
 IMPLICIT NONE
 LOGICAL vmc,dmc

 write(6,*)
 write(6,*)'O---------O'
 write(6,*)'| REBLOCK |'
 write(6,*)'O---------O'
 write(6,*)

! What files are present?
 inquire(file='vmc.hist',exist=vmc)
 inquire(file='dmc.hist',exist=dmc)

 if(.not.(vmc.or.dmc))then
  write(6,*)'Sorry, there are no vmc.hist or dmc.hist files to analyse.'
  stop
 endif ! No hist files found.

! Sort out which file to analyse if more than one possibility exists.
 if(dmc)then
  filename='dmc.hist'
 elseif(vmc)then
  filename='vmc.hist'
 else
  write(6,*)'Bug.'
  stop
 endif
 write(6,*)'Data in '//trim(filename)//' will be analysed.'
 write(6,*)

! Read in data from the file.
 call read_data(dmc)

! Check the data for inconsistencies and get units etc.
 call check_data

! Analyse the data.
 call compute_stats

 write(6,*)'Program finished.'
 write(6,*)

END PROGRAM analyse_qmc
