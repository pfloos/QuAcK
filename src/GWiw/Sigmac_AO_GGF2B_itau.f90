subroutine Sigma_c_GGF2B_he_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B he prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_aecf V_dgbh 
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do fbas=1,nBas2
    do dbas=1,nBas2
     do cbas=1,nBas2
      do ebas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       Binter(abas,ebas,dbas,gbas)=Binter(abas,ebas,dbas,gbas)+G_ao2(fbas,gbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do dbas=1,nBas2
     do ebas=1,nBas2
      do hbas=1,nBas2
       Cinter(abas,hbas,dbas,gbas)=Cinter(abas,hbas,dbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,ebas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do gbas=1,nBas2
     do bbas=1,nBas2
      do dbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,hbas,dbas,gbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c+0.5d0*Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_he_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B he prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = - \sum_cdefgh Gcd Gfh Gge V_aecf V_dgbh 
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do fbas=1,nBas2
    do dbas=1,nBas2
     do cbas=1,nBas2
      do ebas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       Binter(abas,ebas,dbas,hbas)=Binter(abas,ebas,dbas,hbas)+G_ao2(fbas,hbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do dbas=1,nBas2
     do ebas=1,nBas2
      do gbas=1,nBas2
       Cinter(abas,gbas,dbas,hbas)=Cinter(abas,gbas,dbas,hbas)+G_ao3(gbas,ebas)*Binter(abas,ebas,dbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do bbas=1,nBas2
     do dbas=1,nBas2
      do gbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,gbas,dbas,hbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c-Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_eh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B eh prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_cfae V_bhdg
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do cbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do gbas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Binter(dbas,gbas,abas,ebas)=Binter(dbas,gbas,abas,ebas)+G_ao2(fbas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do hbas=1,nBas2
     do gbas=1,nBas2
      do dbas=1,nBas2
       Cinter(dbas,gbas,abas,hbas)=Cinter(dbas,gbas,abas,hbas)+G_ao3(hbas,ebas)*Binter(dbas,gbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do hbas=1,nBas2
     do dbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,gbas,abas,hbas)*db_ERI_AO(bbas,hbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c+0.5d0*Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_eh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B eh prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = - \sum_cdefgh Gcd Geg Ghf V_cfae V_bhdg 
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do cbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do ebas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Binter(dbas,fbas,abas,gbas)=Binter(dbas,fbas,abas,gbas)+G_ao2(ebas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do fbas=1,nBas2
     do hbas=1,nBas2
      do dbas=1,nBas2
       Cinter(dbas,hbas,abas,gbas)=Cinter(dbas,hbas,abas,gbas)+G_ao3(hbas,fbas)*Binter(dbas,fbas,abas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do dbas=1,nBas2
     do hbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,hbas,abas,gbas)*db_ERI_AO(bbas,hbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c-Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_hh_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B hh prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = - \sum_cdefgh Gcd Gfh Gge V_aecf V_hbgd 
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do fbas=1,nBas2
    do dbas=1,nBas2
     do cbas=1,nBas2
      do ebas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       Binter(abas,ebas,dbas,hbas)=Binter(abas,ebas,dbas,hbas)+G_ao2(fbas,hbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do dbas=1,nBas2
     do ebas=1,nBas2
      do gbas=1,nBas2
       Cinter(abas,gbas,dbas,hbas)=Cinter(abas,gbas,dbas,hbas)+G_ao3(gbas,ebas)*Binter(abas,ebas,dbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do dbas=1,nBas2
    do hbas=1,nBas2
     do gbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,gbas,dbas,hbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c-Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_hh_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B hh prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_aecf V_hbgd
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do fbas=1,nBas2
    do dbas=1,nBas2
     do cbas=1,nBas2
      do ebas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       Binter(abas,ebas,dbas,gbas)=Binter(abas,ebas,dbas,gbas)+G_ao2(fbas,gbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do dbas=1,nBas2
     do ebas=1,nBas2
      do hbas=1,nBas2
       Cinter(abas,hbas,dbas,gbas)=Cinter(abas,hbas,dbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,ebas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do dbas=1,nBas2
    do gbas=1,nBas2
     do bbas=1,nBas2
      do hbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,hbas,dbas,gbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c+0.5d0*Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_ee_prime(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B ee prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = - \sum_cdefgh Gcd Geg Ghf V_cfae V_gdhb 
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do cbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do ebas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Binter(dbas,fbas,abas,gbas)=Binter(dbas,fbas,abas,gbas)+G_ao2(ebas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do gbas=1,nBas2
    do fbas=1,nBas2
     do hbas=1,nBas2
      do dbas=1,nBas2
       Cinter(dbas,hbas,abas,gbas)=Cinter(dbas,hbas,abas,gbas)+G_ao3(hbas,fbas)*Binter(dbas,fbas,abas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do bbas=1,nBas2
    do gbas=1,nBas2
     do hbas=1,nBas2
      do dbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,hbas,abas,gbas)*db_ERI_AO(gbas,dbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c-Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

subroutine Sigma_c_GGF2B_ee_prime2(nBas2,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c)

! Gen scGF2B ee prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_ao1(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao2(nBas2,nBas2)
  complex*16,intent(in)         :: G_ao3(nBas2,nBas2)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas
  complex*16,allocatable        :: Sigma_c_tmp(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas2,nBas2)
  complex*16,intent(inout)      :: Ainter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Binter(nBas2,nBas2,nBas2,nBas2)
  complex*16,intent(inout)      :: Cinter(nBas2,nBas2,nBas2,nBas2)

  allocate(Sigma_c_tmp(nBas2,nBas2))
  Sigma_c_tmp=czero

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,ebas,fbas,gbas,hbas)                 & 
  !$OMP          SHARED(nBas2,G_ao1,G_ao2,G_ao3,db_ERI_AO,Sigma_c_tmp,Ainter,Binter,Cinter)

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Geh Ggf V_cfae V_gdhb
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do ebas=1,nBas2
    do fbas=1,nBas2
     do dbas=1,nBas2
      do cbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do ebas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Binter(dbas,fbas,abas,hbas)=Binter(dbas,fbas,abas,hbas)+G_ao2(ebas,hbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do fbas=1,nBas2
     do gbas=1,nBas2
      do dbas=1,nBas2
       Cinter(dbas,gbas,abas,hbas)=Cinter(dbas,gbas,abas,hbas)+G_ao3(gbas,fbas)*Binter(dbas,fbas,abas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do hbas=1,nBas2
    do bbas=1,nBas2
     do gbas=1,nBas2
      do dbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,gbas,abas,hbas)*db_ERI_AO(gbas,dbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  Sigma_c=Sigma_c+0.5d0*Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine

