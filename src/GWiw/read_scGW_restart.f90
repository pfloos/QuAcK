subroutine read_scGW_restart(nBas,nfreqs,ntimes_twice,chem_pot,P_ao,P_ao_hf,G_ao_iw_hf,G_ao_itau,G_ao_itau_hf)

! Read grids for scGW

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)               :: nBas
  integer,intent(in)               :: nfreqs
  integer,intent(in)               :: ntimes_twice
! Local variables
  integer                          :: iunit=312,itau,ifreq,ibas,jbas
  double precision                 :: val_print_r
  complex*16                       :: val_print_c
  
! Output variables
  double precision,intent(out)     :: chem_pot
  double precision,intent(out)     :: P_ao(nBas,nBas)
  double precision,intent(out)     :: P_ao_hf(nBas,nBas)
  complex*16,intent(out)           :: G_ao_iw_hf(nfreqs,nBas,nBas)
  complex*16,intent(out)           :: G_ao_itau(ntimes_twice,nBas,nBas)
  complex*16,intent(out)           :: G_ao_itau_hf(ntimes_twice,nBas,nBas)

  write(*,*)
  write(*,'(a)') ' Reading restart files'
  write(*,*)
  open(unit=iunit,form='unformatted',file='scGW_Gitau_bin',status='old')
  write(*,'(a)') ' Reading scGW_Gitau_bin'
  read(iunit) ibas
  read(iunit) ibas
  do itau=1,ntimes_twice
   do ibas=1,nBas
    do jbas=1,nBas
     read(iunit) val_print_c
     G_ao_itau(itau,ibas,jbas)=val_print_c
    enddo
   enddo
  enddo
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_Pao_bin',status='old')
  write(*,'(a)') ' Reading scGW_Pao_bin'
  read(iunit) ibas
  do ibas=1,nBas
   do jbas=1,nBas
    read(iunit) val_print_r
    P_ao(ibas,jbas)=val_print_r
   enddo
  enddo
  close(iunit)
  open(unit=iunit,form='unformatted',file='scGW_chem_pot_bin',status='old')
  write(*,'(a)') ' Reading scGW_chem_pot_bin'
  read(iunit) chem_pot
  close(iunit)
!  open(unit=iunit,form='unformatted',file='scGW_Gitau_sd_bin',status='old')
!  write(*,'(a)') ' Reading scGW_Gitau_sd_bin'
!  read(iunit) ibas
!  read(iunit) ibas
!  do itau=1,ntimes_twice
!   do ibas=1,nBas
!    do jbas=1,nBas
!     read(iunit) val_print_c
!     G_ao_itau_hf(itau,ibas,jbas)=val_print_c
!    enddo
!   enddo
!  enddo
!  close(iunit)
!  open(unit=iunit,form='unformatted',file='scGW_Giw_sd_bin',status='old')
!  write(*,'(a)') ' Reading scGW_Giw_sd_bin'
!  read(iunit) ibas
!  read(iunit) ibas
!  do ifreq=1,nfreqs
!   do ibas=1,nBas
!    do jbas=1,nBas
!     read(iunit) val_print_c
!     G_ao_iw_hf(ifreq,ibas,jbas)=val_print_c
!    enddo
!   enddo
!  enddo
!  close(iunit)
!  open(unit=iunit,form='unformatted',file='scGW_Pao_sd_bin',status='old')
!  write(*,'(a)') ' Reading scGW_Pao_sd_bin'
!  read(iunit) ibas
!  do ibas=1,nBas
!   do jbas=1,nBas
!    read(iunit) val_print_r
!    P_ao_hf(ibas,jbas)=val_print_r
!   enddo
!  enddo
!  close(iunit)

end subroutine
