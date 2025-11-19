subroutine write_scGWB_restart(nBas_twice,ntimes,ntimes_twice,nfreqs,chem_pot,R_ao,R_ao_hfb,G_ao_itau,G_ao_itau_hfb, &
                               G_ao_iw_hfb,DeltaG_ao_iw)

! Read grids for scGWB

  implicit none
  include 'parameters.h'

! Input variables
 integer,intent(in)              :: nBas_twice
 integer,intent(in)              :: ntimes
 integer,intent(in)              :: nfreqs
 integer,intent(in)              :: ntimes_twice
 double precision,intent(in)     :: chem_pot
 double precision,intent(in)     :: R_ao(nBas_twice,nBas_twice)
 double precision,intent(in)     :: R_ao_hfb(nBas_twice,nBas_twice)
 complex*16,intent(in)           :: G_ao_iw_hfb(nfreqs,nBas_twice,nBas_twice)
 complex*16,intent(in)           :: DeltaG_ao_iw(nfreqs,nBas_twice,nBas_twice)
 complex*16,intent(in)           :: G_ao_itau(ntimes_twice,nBas_twice,nBas_twice)
 complex*16,intent(in)           :: G_ao_itau_hfb(ntimes_twice,nBas_twice,nBas_twice)

! Local variables
 integer                         :: iunit=312,itau,ifreq,ibas,jbas
 double precision                :: val_print_r
 complex*16                      :: val_print_c

  open(unit=iunit,form='unformatted',file='scGWB_Gitau_bin')
  write(iunit) nBas_twice
  write(iunit) ntimes
  do itau=1,ntimes_twice
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     val_print_c = G_ao_itau(itau,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Giw_bin')
  write(iunit) nBas_twice
  write(iunit) ntimes
  do ifreq=1,nfreqs
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     val_print_c = G_ao_iw_hfb(ifreq,ibas,jbas)+DeltaG_ao_iw(ifreq,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Rao_bin')
  write(iunit) nBas_twice
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    val_print_r = R_ao(ibas,jbas)
    if(abs(val_print_r)<1d-8) val_print_r=czero
    write(iunit) val_print_r
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_chem_pot_bin')
  write(iunit) chem_pot
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Gitau_hfb_bin')
  write(iunit) nBas_twice
  write(iunit) ntimes
  do itau=1,ntimes_twice
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     val_print_c = G_ao_itau_hfb(itau,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Giw_hfb_bin')
  write(iunit) nBas_twice
  write(iunit) ntimes
  do ifreq=1,nfreqs
   do ibas=1,nBas_twice
    do jbas=1,nBas_twice
     val_print_c = G_ao_iw_hfb(ifreq,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGWB_Rao_hfb_bin')
  write(iunit) nBas_twice
  do ibas=1,nBas_twice
   do jbas=1,nBas_twice
    val_print_r = R_ao_hfb(ibas,jbas)
    if(abs(val_print_r)<1d-8) val_print_r=czero
    write(iunit) val_print_r
   enddo
  enddo
  write(iunit) iunit
  close(iunit)

end subroutine
