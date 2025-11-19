subroutine read_scGWB_restart(nBas_twice,nfreqs,ntimes_twice,chem_pot,R_ao,R_ao_hfb,G_ao_iw_hfb,G_ao_itau,G_ao_itau_hfb,read_SD_chkp)

! Read grids for scGWB

  implicit none
  include 'parameters.h'

! Input variables
  logical,intent(in)               :: read_SD_chkp
  integer,intent(in)               :: nBas_twice
  integer,intent(in)               :: nfreqs
  integer,intent(in)               :: ntimes_twice
! Local variables
  integer                          :: iunit=312,itau,ifreq,ibas,jbas
  double precision                 :: val_print_r
  complex*16                       :: val_print_c
  
! Output variables
  double precision,intent(out)     :: chem_pot
  double precision,intent(out)     :: R_ao(nBas_twice,nBas_twice)
  double precision,intent(out)     :: R_ao_hfb(nBas_twice,nBas_twice)
  complex*16,intent(out)           :: G_ao_iw_hfb(nfreqs,nBas_twice,nBas_twice)
  complex*16,intent(out)           :: G_ao_itau(ntimes_twice,nBas_twice,nBas_twice)
  complex*16,intent(out)           :: G_ao_itau_hfb(ntimes_twice,nBas_twice,nBas_twice)

  write(*,*)
  write(*,'(a)') ' Reading restart files'
  write(*,*)
  open(unit=iunit,form='unformatted',file='scGWB_Gitau_bin',status='old')
  write(*,'(a)') ' Reading scGWB_Gitau_bin'
  read(iunit) ibas
  read(iunit) ibas
  do itau=1,ntimes_twice
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     read(iunit) val_print_c
     G_ao_itau(itau,ibas,jbas)=val_print_c
    enddo
   enddo
  enddo
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Rao_bin',status='old')
  write(*,'(a)') ' Reading scGWB_Rao_bin'
  read(iunit) ibas
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    read(iunit) val_print_r
    R_ao(ibas,jbas)=val_print_r
   enddo
  enddo
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_chem_pot_bin',status='old')
  write(*,'(a)') ' Reading scGWB_chem_pot_bin'
  read(iunit) chem_pot
  close(iunit)
  if(read_SD_chkp) then
   open(unit=iunit,form='unformatted',file='scGWB_Gitau_hfb_bin',status='old')
   write(*,'(a)') ' Reading scGWB_Gitau_hfb_bin'
   read(iunit) ibas
   read(iunit) ibas
   do itau=1,ntimes_twice
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      read(iunit) val_print_c
      G_ao_itau_hfb(itau,ibas,jbas)=val_print_c
     enddo
    enddo
   enddo
   close(iunit)
   open(unit=iunit,form='unformatted',file='scGWB_Giw_hfb_bin',status='old')
   write(*,'(a)') ' Reading scGWB_Giw_hfb_bin'
   read(iunit) ibas
   read(iunit) ibas
   do ifreq=1,nfreqs
    do ibas=1,nBas_twice
     do jbas=1,nBas_twice
      read(iunit) val_print_c
      G_ao_iw_hfb(ifreq,ibas,jbas)=val_print_c
     enddo
    enddo
   enddo
   close(iunit)
   open(unit=iunit,form='unformatted',file='scGWB_Rao_hfb_bin',status='old')
   write(*,'(a)') ' Reading scGWB_Rao_hfb_bin'
   read(iunit) ibas
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     read(iunit) val_print_r
     R_ao_hfb(ibas,jbas)=val_print_r
    enddo
   enddo
   close(iunit)
  endif

end subroutine
