subroutine read_scGW_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                           cosw2t_weight,sinw2t_weight,verbose)

! Read grids for scGW

  implicit none
  include 'parameters.h'

! Input variables
 integer,intent(in)              :: verbose
 integer,intent(in)              :: ntimes
 integer,intent(in)              :: nfreqs

! Local variables
 logical                         :: file_exists
 integer                         :: iunit=312
 integer                         :: itau,ifreq,ivar
 double precision                :: error_I
 double precision,allocatable    :: I_weight(:,:)

! Output variables
 double precision,intent(out)    :: tcoord(ntimes)
 double precision,intent(out)    :: tweight(ntimes)
 double precision,intent(out)    :: wcoord(nfreqs)
 double precision,intent(out)    :: wweight(nfreqs)
 double precision,intent(out)    :: sint2w_weight(nfreqs,ntimes)
 double precision,intent(out)    :: cost2w_weight(nfreqs,ntimes)
 double precision,intent(out)    :: cosw2t_weight(ntimes,nfreqs)
 double precision,intent(out)    :: sinw2t_weight(ntimes,nfreqs)

 allocate(I_weight(nfreqs,ntimes))
 if(verbose/=0)  write(*,*)
 inquire(file='./grids/tcoord.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading tcoord from tcoord.txt'
  tcoord=0d0
  open(unit=iunit, form='formatted', file='./grids/tcoord.txt', status='old')
  read(iunit,*) ivar
  do itau=1,ntimes
   read(iunit,*) tcoord(itau)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the tcoord.txt file'
  return
 endif
 inquire(file='./grids/tweight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading tweight from tweight.txt'
  tweight=0d0
  open(unit=iunit, form='formatted', file='./grids/tweight.txt', status='old')
  read(iunit,*) ivar
  do itau=1,ntimes
   read(iunit,*) tweight(itau)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the tweight.txt file'
  return
 endif
 inquire(file='./grids/wcoord.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading wcoord from wcoord.txt'
  wcoord=0d0
  open(unit=iunit, form='formatted', file='./grids/wcoord.txt', status='old')
  read(iunit,*) ivar
  do ifreq=1,nfreqs
   read(iunit,*) wcoord(ifreq)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the wcoord.txt file'
  return
 endif
 inquire(file='./grids/wweight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading wweight from wweight.txt'
  wweight=0d0
  open(unit=iunit, form='formatted', file='./grids/wweight.txt', status='old')
  read(iunit,*) ivar
  do ifreq=1,nfreqs
   read(iunit,*) wweight(ifreq)
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the wweight.txt file'
  return
 endif
 inquire(file='./grids/cost2w_weight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading cost2w_weight from cost2w_weight.txt'
  cost2w_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/cost2w_weight.txt', status='old')
  read(iunit,*) ivar
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) cost2w_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the cost2w_weight.txt file'
  return
 endif
 inquire(file='./grids/cosw2t_weight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading cosw2t_weight from cosw2t_weight.txt'
  cosw2t_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/cosw2t_weight.txt', status='old')
  read(iunit,*) ivar
  do itau=1,ntimes
   do ifreq=1,nfreqs
    read(iunit,*) cosw2t_weight(itau,ifreq)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the cosw2t_weight.txt file'
  return
 endif
 inquire(file='./grids/sint2w_weight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading sint2w_weight from sint2w_weight.txt'
  sint2w_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/sint2w_weight.txt', status='old')
  read(iunit,*) ivar
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) sint2w_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the sint2w_weight.txt file'
  return
 endif
 inquire(file='./grids/sinw2t_weight.txt', exist=file_exists)
 if(file_exists) then
  if(verbose/=0) write(*,*) 'Reading sinw2t_weight from sinw2t_weight.txt'
  sinw2t_weight=0d0
  open(unit=iunit, form='formatted', file='./grids/sinw2t_weight.txt', status='old')
  read(iunit,*) ivar
  do ifreq=1,nfreqs
   do itau=1,ntimes
    read(iunit,*) sinw2t_weight(ifreq,itau)
   enddo
  enddo
  close(iunit)
 else
  write(*,*) ' Error! Could not find the sinw2t_weight.txt file'
  return
 endif
 if(verbose/=0)  write(*,'(a,i5,a)') ' Using ',nfreqs,' frequencies and times'
 I_weight=matmul(cost2w_weight,cosw2t_weight)
 do ifreq=1,nfreqs
  I_weight(ifreq,ifreq)=I_weight(ifreq,ifreq)-1d0
 enddo
 error_I=0d0
 do ifreq=1,nfreqs
  do itau=1,ntimes
   error_I=error_I+abs(I_weight(ifreq,itau))
  enddo
 enddo
 if(verbose/=0)  write(*,'(a,f20.5)') ' Deviation from Identity in Cos-Cos ',error_I
 I_weight=matmul(sint2w_weight,sinw2t_weight)
 do ifreq=1,nfreqs
  I_weight(ifreq,ifreq)=I_weight(ifreq,ifreq)-1d0
 enddo
 error_I=0d0
 do ifreq=1,nfreqs
  do itau=1,ntimes
   error_I=error_I+abs(I_weight(ifreq,itau))
  enddo
 enddo
 if(verbose/=0)  write(*,'(a,f20.5)') ' Deviation from Identity in Sin-Sin ',error_I
 if(verbose/=0)  write(*,*)
 deallocate(I_weight)

end subroutine

