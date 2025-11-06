subroutine write_scGW_restart(nBas,ntimes,ntimes_twice,nfreqs,chem_pot,P_ao,P_ao_hf,G_ao_itau,G_ao_itau_hf, &
                              G_ao_iw_hf,DeltaG_ao_iw)

! Read grids for scGW

  implicit none
  include 'parameters.h'

! Input variables
 integer,intent(in)              :: nBas
 integer,intent(in)              :: ntimes
 integer,intent(in)              :: nfreqs
 integer,intent(in)              :: ntimes_twice
 double precision,intent(in)     :: chem_pot
 double precision,intent(in)     :: P_ao(nBas,nBas)
 double precision,intent(in)     :: P_ao_hf(nBas,nBas)
 complex*16,intent(in)           :: G_ao_iw_hf(nfreqs,nBas,nBas)
 complex*16,intent(in)           :: DeltaG_ao_iw(nfreqs,nBas,nBas)
 complex*16,intent(in)           :: G_ao_itau(ntimes_twice,nBas,nBas)
 complex*16,intent(in)           :: G_ao_itau_hf(ntimes_twice,nBas,nBas)

! Local variables
 integer                         :: iunit=312,itau,ifreq,ibas,jbas
 double precision                :: val_print_r
 complex*16                      :: val_print_c

  open(unit=iunit,form='unformatted',file='scGW_Gitau_bin')
  write(iunit) nBas
  write(iunit) ntimes
  do itau=1,ntimes_twice
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_itau(itau,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Giw_bin')
  write(iunit) nBas
  write(iunit) ntimes
  do ifreq=1,nfreqs
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_iw_hf(ifreq,ibas,jbas)+DeltaG_ao_iw(ifreq,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Pao_bin')
  write(iunit) nBas
  do ibas=1,nBas
   do jbas=1,nBas
    val_print_r = P_ao(ibas,jbas)
    if(abs(val_print_r)<1d-8) val_print_r=czero
    write(iunit) val_print_r
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_chem_pot_bin')
  write(iunit) chem_pot
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Gitau_sd_bin')
  write(iunit) nBas
  write(iunit) ntimes
  do itau=1,ntimes_twice
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_itau_hf(itau,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Giw_sd_bin')
  write(iunit) nBas
  write(iunit) ntimes
  do ifreq=1,nfreqs
   do ibas=1,nBas
    do jbas=1,nBas
     val_print_c = G_ao_iw_hf(ifreq,ibas,jbas)
     if(abs(val_print_c)<1d-8) val_print_c=czero
     write(iunit) val_print_c
    enddo
   enddo
  enddo
  write(iunit) iunit
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Pao_sd_bin')
  write(iunit) nBas
  do ibas=1,nBas
   do jbas=1,nBas
    val_print_r = P_ao_hf(ibas,jbas)
    if(abs(val_print_r)<1d-8) val_print_r=czero
    write(iunit) val_print_r
   enddo
  enddo
  write(iunit) iunit
  close(iunit)

end subroutine
