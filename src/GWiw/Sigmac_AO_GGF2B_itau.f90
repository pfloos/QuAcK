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

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_aecf V_dgbh 
  do abas=1,nBas2
   do ebas=1,nBas2
    do cbas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do fbas=1,nBas2
      do gbas=1,nBas2
       Binter(abas,ebas,dbas,gbas)=Binter(abas,ebas,dbas,gbas)+G_ao2(fbas,gbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do gbas=1,nBas2
      do hbas=1,nBas2
       Cinter(abas,hbas,dbas,gbas)=Cinter(abas,hbas,dbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,ebas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do hbas=1,nBas2
    do dbas=1,nBas2
     do gbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,hbas,dbas,gbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = - \sum_cdefgh Gcd Gfh Gge V_aecf V_dgbh 
  do abas=1,nBas2
   do ebas=1,nBas2
    do cbas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do fbas=1,nBas2
      do hbas=1,nBas2
       Binter(abas,ebas,dbas,hbas)=Binter(abas,ebas,dbas,hbas)+G_ao2(fbas,hbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do hbas=1,nBas2
      do gbas=1,nBas2
       Cinter(abas,gbas,dbas,hbas)=Cinter(abas,gbas,dbas,hbas)+G_ao3(gbas,ebas)*Binter(abas,ebas,dbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do gbas=1,nBas2
    do dbas=1,nBas2
     do hbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,gbas,dbas,hbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_cfae V_bhdg
  do cbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do dbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do gbas=1,nBas2
       Binter(dbas,gbas,abas,ebas)=Binter(dbas,gbas,abas,ebas)+G_ao2(fbas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do gbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do hbas=1,nBas2
       Cinter(dbas,gbas,abas,hbas)=Cinter(dbas,gbas,abas,hbas)+G_ao3(hbas,ebas)*Binter(dbas,gbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do gbas=1,nBas2
    do abas=1,nBas2
     do hbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,gbas,abas,hbas)*db_ERI_AO(bbas,hbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = - \sum_cdefgh Gcd Geg Ghf V_cfae V_bhdg 
  do cbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do dbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do gbas=1,nBas2
       Binter(dbas,fbas,abas,gbas)=Binter(dbas,fbas,abas,gbas)+G_ao2(ebas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do gbas=1,nBas2
      do hbas=1,nBas2
       Cinter(dbas,hbas,abas,gbas)=Cinter(dbas,hbas,abas,gbas)+G_ao3(hbas,fbas)*Binter(dbas,fbas,abas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do hbas=1,nBas2
    do abas=1,nBas2
     do gbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,hbas,abas,gbas)*db_ERI_AO(bbas,hbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = - \sum_cdefgh Gcd Gfh Gge V_aecf V_hbgd 
  do abas=1,nBas2
   do ebas=1,nBas2
    do cbas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do fbas=1,nBas2
      do hbas=1,nBas2
       Binter(abas,ebas,dbas,hbas)=Binter(abas,ebas,dbas,hbas)+G_ao2(fbas,hbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do hbas=1,nBas2
      do gbas=1,nBas2
       Cinter(abas,gbas,dbas,hbas)=Cinter(abas,gbas,dbas,hbas)+G_ao3(gbas,ebas)*Binter(abas,ebas,dbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do gbas=1,nBas2
    do dbas=1,nBas2
     do hbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,gbas,dbas,hbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Gfg Ghe V_aecf V_hbgd
  do abas=1,nBas2
   do ebas=1,nBas2
    do cbas=1,nBas2
     do fbas=1,nBas2
      do dbas=1,nBas2
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*db_ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do fbas=1,nBas2
      do gbas=1,nBas2
       Binter(abas,ebas,dbas,gbas)=Binter(abas,ebas,dbas,gbas)+G_ao2(fbas,gbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do ebas=1,nBas2
    do dbas=1,nBas2
     do gbas=1,nBas2
      do hbas=1,nBas2
       Cinter(abas,hbas,dbas,gbas)=Cinter(abas,hbas,dbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,ebas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas2
   do hbas=1,nBas2
    do dbas=1,nBas2
     do gbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(abas,hbas,dbas,gbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = - \sum_cdefgh Gcd Geg Ghf V_cfae V_gdhb 
  do cbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do dbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do gbas=1,nBas2
       Binter(dbas,fbas,abas,gbas)=Binter(dbas,fbas,abas,gbas)+G_ao2(ebas,gbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do gbas=1,nBas2
      do hbas=1,nBas2
       Cinter(dbas,hbas,abas,gbas)=Cinter(dbas,hbas,abas,gbas)+G_ao3(hbas,fbas)*Binter(dbas,fbas,abas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do hbas=1,nBas2
    do abas=1,nBas2
     do gbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,hbas,abas,gbas)*db_ERI_AO(gbas,dbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo

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

  ! Sigma_c_ab = 0.5 \sum_cdefgh Gcd Geh Ggf V_cfae V_gdhb
  do cbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do dbas=1,nBas2
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*db_ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do ebas=1,nBas2
      do hbas=1,nBas2
       Binter(dbas,fbas,abas,hbas)=Binter(dbas,fbas,abas,hbas)+G_ao2(ebas,hbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do fbas=1,nBas2
    do abas=1,nBas2
     do hbas=1,nBas2
      do gbas=1,nBas2
       Cinter(dbas,gbas,abas,hbas)=Cinter(dbas,gbas,abas,hbas)+G_ao3(gbas,fbas)*Binter(dbas,fbas,abas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas2
   do gbas=1,nBas2
    do abas=1,nBas2
     do hbas=1,nBas2
      do bbas=1,nBas2
       Sigma_c_tmp(abas,bbas)=Sigma_c_tmp(abas,bbas)+Cinter(dbas,gbas,abas,hbas)*db_ERI_AO(gbas,dbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo

  Sigma_c=Sigma_c+0.5d0*Sigma_c_tmp

  deallocate(Sigma_c_tmp)

end subroutine


subroutine Sigma_c_GGF2B_brut(nBas2,nBas4,G_plus,G_minus,db_ERI_AO,Sigma_c_plus,Sigma_c_minus) 

! Restricted scGF2B all block M^8

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nBas4

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_plus(nBas4,nBas4)
  complex*16,intent(in)         :: G_minus(nBas4,nBas4)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

  complex*16                    :: Integ_val
  complex*16,allocatable        :: G_he_plus(:,:)
  complex*16,allocatable        :: G_hh_plus(:,:)
  complex*16,allocatable        :: G_ee_plus(:,:)
  complex*16,allocatable        :: G_eh_plus(:,:)
  complex*16,allocatable        :: G_he_minus(:,:)
  complex*16,allocatable        :: G_hh_minus(:,:)
  complex*16,allocatable        :: G_ee_minus(:,:)
  complex*16,allocatable        :: G_eh_minus(:,:)
  complex*16,allocatable        :: Sigma_he_plus(:,:)
  complex*16,allocatable        :: Sigma_hh_plus(:,:)
  complex*16,allocatable        :: Sigma_ee_plus(:,:)
  complex*16,allocatable        :: Sigma_eh_plus(:,:)
  complex*16,allocatable        :: Sigma_he_minus(:,:)
  complex*16,allocatable        :: Sigma_hh_minus(:,:)
  complex*16,allocatable        :: Sigma_ee_minus(:,:)
  complex*16,allocatable        :: Sigma_eh_minus(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c_plus(nBas4,nBas4)
  complex*16,intent(inout)      :: Sigma_c_minus(nBas4,nBas4)

  ! Allocate arrays

  allocate(G_he_plus(nBas2,nBas2)) 
  allocate(G_hh_plus(nBas2,nBas2))
  allocate(G_ee_plus(nBas2,nBas2))
  allocate(G_eh_plus(nBas2,nBas2))
  allocate(G_he_minus(nBas2,nBas2))
  allocate(G_hh_minus(nBas2,nBas2))
  allocate(G_ee_minus(nBas2,nBas2))
  allocate(G_eh_minus(nBas2,nBas2))
  allocate(Sigma_he_plus(nBas2,nBas2)) 
  allocate(Sigma_hh_plus(nBas2,nBas2))
  allocate(Sigma_ee_plus(nBas2,nBas2))
  allocate(Sigma_eh_plus(nBas2,nBas2))
  allocate(Sigma_he_minus(nBas2,nBas2))
  allocate(Sigma_hh_minus(nBas2,nBas2))
  allocate(Sigma_ee_minus(nBas2,nBas2))
  allocate(Sigma_eh_minus(nBas2,nBas2))

  ! Initialize 
  
  G_he_plus(1:nBas2,1:nBas2)=G_plus(      1:nBas2,      1:nBas2)
  G_hh_plus(1:nBas2,1:nBas2)=G_plus(      1:nBas2,nBas2+1:nBas4)    
  G_ee_plus(1:nBas2,1:nBas2)=G_plus(nBas2+1:nBas4,      1:nBas2)   
  G_eh_plus(1:nBas2,1:nBas2)=G_plus(nBas2+1:nBas4,nBas2+1:nBas4)  
  G_he_minus(1:nBas2,1:nBas2)=G_minus(      1:nBas2,      1:nBas2)
  G_hh_minus(1:nBas2,1:nBas2)=G_minus(      1:nBas2,nBas2+1:nBas4)    
  G_ee_minus(1:nBas2,1:nBas2)=G_minus(nBas2+1:nBas4,      1:nBas2)   
  G_eh_minus(1:nBas2,1:nBas2)=G_minus(nBas2+1:nBas4,nBas2+1:nBas4)  
  Sigma_he_plus=czero
  Sigma_hh_plus=czero
  Sigma_ee_plus=czero
  Sigma_eh_plus=czero
  Sigma_he_minus=czero
  Sigma_hh_minus=czero
  Sigma_ee_minus=czero
  Sigma_eh_minus=czero

  ! M^8 loop
  !write(*,*)
  !write(*,*) 'Computing Bog. Sigma_c (M^8)'
  !write(*,*)

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas,Integ_val)       &
  !$OMP          SHARED(nBas2,nBas4,G_plus,G_minus,db_ERI_AO,                    &
  !$OMP &        G_he_plus,G_hh_plus,G_ee_plus,G_eh_plus,G_he_minus,G_hh_minus,   &
  !$OMP &        G_ee_minus,G_eh_minus,Sigma_he_plus,Sigma_hh_plus,Sigma_ee_plus, &
  !$OMP &        Sigma_eh_plus,Sigma_he_minus,Sigma_hh_minus,Sigma_ee_minus,Sigma_eh_minus)
  block ! Use block to define local arrays
  complex*16,allocatable        :: Sigma_he_plus_th(:,:)
  complex*16,allocatable        :: Sigma_hh_plus_th(:,:)
  complex*16,allocatable        :: Sigma_ee_plus_th(:,:)
  complex*16,allocatable        :: Sigma_eh_plus_th(:,:)
  complex*16,allocatable        :: Sigma_he_minus_th(:,:)
  complex*16,allocatable        :: Sigma_hh_minus_th(:,:)
  complex*16,allocatable        :: Sigma_ee_minus_th(:,:)
  complex*16,allocatable        :: Sigma_eh_minus_th(:,:)

  ! Allocate local arrays
  allocate(Sigma_he_plus_th(nBas2,nBas2)) 
  allocate(Sigma_hh_plus_th(nBas2,nBas2))
  allocate(Sigma_ee_plus_th(nBas2,nBas2))
  allocate(Sigma_eh_plus_th(nBas2,nBas2))
  allocate(Sigma_he_minus_th(nBas2,nBas2))
  allocate(Sigma_hh_minus_th(nBas2,nBas2))
  allocate(Sigma_ee_minus_th(nBas2,nBas2))
  allocate(Sigma_eh_minus_th(nBas2,nBas2))

  ! Initialize local arrays
  Sigma_he_plus_th=czero
  Sigma_hh_plus_th=czero
  Sigma_ee_plus_th=czero
  Sigma_eh_plus_th=czero
  Sigma_he_minus_th=czero
  Sigma_hh_minus_th=czero
  Sigma_ee_minus_th=czero
  Sigma_eh_minus_th=czero

  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do bbas=1,nBas2
    do cbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       do fbas=1,nBas2
        do gbas=1,nBas2
         do hbas=1,nBas2
           ! he 2'
           Integ_val=0.5d0*db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
           Sigma_he_plus_th(abas,bbas) = Sigma_he_plus_th(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_he_plus(fbas,gbas) *G_he_minus(hbas,ebas)            
           Sigma_he_minus_th(abas,bbas)=Sigma_he_minus_th(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_he_minus(fbas,gbas)*G_he_plus(hbas,ebas)            
           ! he 2''
           Integ_val=-db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
           Sigma_he_plus_th(abas,bbas) = Sigma_he_plus_th(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_hh_plus(fbas,hbas) *G_ee_minus(gbas,ebas)            
           Sigma_he_minus_th(abas,bbas)=Sigma_he_minus_th(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_hh_minus(fbas,hbas)*G_ee_plus(gbas,ebas) 
           ! eh 2'
           Integ_val=0.5d0*db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(bbas,hbas,dbas,gbas)
           Sigma_eh_plus_th(abas,bbas) = Sigma_eh_plus_th(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_eh_plus(fbas,gbas) *G_eh_minus(hbas,ebas)            
           Sigma_eh_minus_th(abas,bbas)=Sigma_eh_minus_th(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_eh_minus(fbas,gbas)*G_eh_plus(hbas,ebas)            
           ! eh 2''
           Integ_val=-db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(bbas,hbas,dbas,gbas)
           Sigma_eh_plus_th(abas,bbas) = Sigma_eh_plus_th(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_hh_plus(ebas,gbas) *G_ee_minus(hbas,fbas)            
           Sigma_eh_minus_th(abas,bbas)=Sigma_eh_minus_th(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_hh_minus(ebas,gbas)*G_ee_plus(hbas,fbas) 
           ! hh 2'
           Integ_val=-db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
           Sigma_hh_plus_th(abas,bbas) = Sigma_hh_plus_th(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_he_plus(fbas,hbas) *G_he_minus(gbas,ebas)            
           Sigma_hh_minus_th(abas,bbas)=Sigma_hh_minus_th(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_he_minus(fbas,hbas)*G_he_plus(gbas,ebas) 
           ! hh 2''
           Integ_val=0.5d0*db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
           Sigma_hh_plus_th(abas,bbas) = Sigma_hh_plus_th(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_hh_plus(fbas,gbas) *G_ee_minus(hbas,ebas)            
           Sigma_hh_minus_th(abas,bbas)=Sigma_hh_minus_th(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_hh_minus(fbas,gbas)*G_ee_plus(hbas,ebas)            
           ! ee 2'
           Integ_val=-db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(gbas,dbas,hbas,bbas)
           Sigma_ee_plus_th(abas,bbas) = Sigma_ee_plus_th(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_he_plus(ebas,gbas) *G_he_minus(hbas,fbas)            
           Sigma_ee_minus_th(abas,bbas)=Sigma_ee_minus_th(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_he_minus(ebas,gbas)*G_he_plus(hbas,fbas) 
           ! ee 2''
           Integ_val=0.5*db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(gbas,dbas,hbas,bbas)
           Sigma_ee_plus_th(abas,bbas) = Sigma_ee_plus_th(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_hh_plus(ebas,hbas) *G_ee_minus(gbas,fbas)            
           Sigma_ee_minus_th(abas,bbas)=Sigma_ee_minus_th(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_hh_minus(ebas,hbas)*G_ee_plus(gbas,fbas)            
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO

  ! Add contributions (one thread at a time) 
  !$OMP CRITICAL

  Sigma_he_plus = Sigma_he_plus + Sigma_he_plus_th    
  Sigma_hh_plus = Sigma_hh_plus + Sigma_hh_plus_th
  Sigma_ee_plus = Sigma_ee_plus + Sigma_ee_plus_th
  Sigma_eh_plus = Sigma_eh_plus + Sigma_eh_plus_th
  Sigma_he_minus = Sigma_he_minus + Sigma_he_minus_th
  Sigma_hh_minus = Sigma_hh_minus + Sigma_hh_minus_th
  Sigma_ee_minus = Sigma_ee_minus + Sigma_ee_minus_th
  Sigma_eh_minus = Sigma_eh_minus + Sigma_eh_minus_th

  !$OMP END CRITICAL

  ! Wait for all to finish
  !$OMP BARRIER 

  ! Deallocate local arrays 
  deallocate(Sigma_he_plus_th) 
  deallocate(Sigma_hh_plus_th)
  deallocate(Sigma_ee_plus_th)
  deallocate(Sigma_eh_plus_th)
  deallocate(Sigma_he_minus_th)
  deallocate(Sigma_hh_minus_th)
  deallocate(Sigma_ee_minus_th)
  deallocate(Sigma_eh_minus_th)

  end block !
  !$OMP END PARALLEL

  ! Set Sigma_c Gorkov

  Sigma_c_plus(      1:nBas2,       1:nBas2)=Sigma_he_plus(1:nBas2,1:nBas2)
  Sigma_c_plus(      1:nBas2, nBas2+1:nBas4)=Sigma_hh_plus(1:nBas2,1:nBas2)    
  Sigma_c_plus(nBas2+1:nBas4,       1:nBas2)=Sigma_ee_plus(1:nBas2,1:nBas2)   
  Sigma_c_plus(nBas2+1:nBas4, nBas2+1:nBas4)=Sigma_eh_plus(1:nBas2,1:nBas2)  
  Sigma_c_minus(     1:nBas2,       1:nBas2)=Sigma_he_minus(1:nBas2,1:nBas2)
  Sigma_c_minus(     1:nBas2, nBas2+1:nBas4)=Sigma_hh_minus(1:nBas2,1:nBas2)    
  Sigma_c_minus(nBas2+1:nBas4,      1:nBas2)=Sigma_ee_minus(1:nBas2,1:nBas2)   
  Sigma_c_minus(nBas2+1:nBas4,nBas2+1:nBas4)=Sigma_eh_minus(1:nBas2,1:nBas2)  

  ! Deallocate arrays 

  deallocate(G_he_plus)
  deallocate(G_hh_plus)
  deallocate(G_ee_plus)
  deallocate(G_eh_plus)
  deallocate(G_he_minus)
  deallocate(G_hh_minus)
  deallocate(G_ee_minus)
  deallocate(G_eh_minus)
  deallocate(Sigma_he_plus) 
  deallocate(Sigma_hh_plus)
  deallocate(Sigma_ee_plus)
  deallocate(Sigma_eh_plus)
  deallocate(Sigma_he_minus)
  deallocate(Sigma_hh_minus)
  deallocate(Sigma_ee_minus)
  deallocate(Sigma_eh_minus)

end subroutine

