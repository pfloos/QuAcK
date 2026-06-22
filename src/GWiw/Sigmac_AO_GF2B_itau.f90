subroutine Sigma_c_GF2B_he_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c) 

! Restricted scGF2B he prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfg Ghe v_aecf (2 v_dgbh - v_dghb )
  do abas=1,nBas
   do cbas=1,nBas
    do ebas=1,nBas
     do gbas=1,nBas
      do fbas=1,nBas
       Ainter(abas,cbas,ebas,gbas)=Ainter(abas,cbas,ebas,gbas)+G_ao2(fbas,gbas)*ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do dbas=1,nBas
    do ebas=1,nBas
     do gbas=1,nBas
      do cbas=1,nBas
       Binter(abas,dbas,ebas,gbas)=Binter(abas,dbas,ebas,gbas)+G_ao1(cbas,dbas)*Ainter(abas,cbas,ebas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do dbas=1,nBas
    do hbas=1,nBas
     do gbas=1,nBas
      do ebas=1,nBas
       Cinter(abas,dbas,hbas,gbas)=Cinter(abas,dbas,hbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,dbas,ebas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do dbas=1,nBas
     do hbas=1,nBas
      do gbas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(abas,dbas,hbas,gbas)*(2d0*ERI_AO(dbas,gbas,bbas,hbas)-ERI_AO(dbas,gbas,hbas,bbas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine 

subroutine Sigma_c_GF2B_he_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c) 

! Restricted scGF2B he prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfh Gge v_dgbh (-2 v_aecf + v_aefc )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(cbas,gbas,bbas,hbas)=Ainter(cbas,gbas,bbas,hbas)+G_ao1(cbas,dbas)*ERI_AO(dbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(cbas,gbas,bbas,fbas)=Binter(cbas,gbas,bbas,fbas)+G_ao2(fbas,hbas)*Ainter(cbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do fbas=1,nBas
       Cinter(cbas,ebas,bbas,fbas)=Cinter(cbas,ebas,bbas,fbas)+G_ao3(gbas,ebas)*Binter(cbas,gbas,bbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do cbas=1,nBas
     do ebas=1,nBas
      do fbas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(cbas,ebas,bbas,fbas)*(-2d0*ERI_AO(abas,ebas,cbas,fbas)+ERI_AO(abas,ebas,fbas,cbas))
      enddo
     enddo
    enddo
   enddo
  enddo

  Ainter=czero; Binter=czero; Cinter=czero;

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfh Gge v_dghb ( v_aecf - v_aefc )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(cbas,gbas,bbas,hbas)=Ainter(cbas,gbas,bbas,hbas)+G_ao1(cbas,dbas)*ERI_AO(dbas,gbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(cbas,gbas,bbas,fbas)=Binter(cbas,gbas,bbas,fbas)+G_ao2(fbas,hbas)*Ainter(cbas,gbas,bbas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do fbas=1,nBas
       Cinter(cbas,ebas,bbas,fbas)=Cinter(cbas,ebas,bbas,fbas)+G_ao3(gbas,ebas)*Binter(cbas,gbas,bbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do cbas=1,nBas
     do ebas=1,nBas
      do fbas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(cbas,ebas,bbas,fbas)*(ERI_AO(abas,ebas,cbas,fbas)-ERI_AO(abas,ebas,fbas,cbas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine 

subroutine Sigma_c_GF2B_eh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c) 

! Restricted scGF2B eh prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Geg Ghf v_bhdg (-2 v_cfae + v_fcae )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(bbas,hbas,cbas,gbas)=Ainter(bbas,hbas,cbas,gbas)+G_ao1(cbas,dbas)*ERI_AO(bbas,hbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(bbas,hbas,cbas,ebas)=Binter(bbas,hbas,cbas,ebas)+G_ao2(ebas,gbas)*Ainter(bbas,hbas,cbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do bbas=1,nBas
    do hbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(bbas,fbas,cbas,ebas)=Cinter(bbas,fbas,cbas,ebas)+G_ao3(hbas,fbas)*Binter(bbas,hbas,cbas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(bbas,fbas,cbas,ebas)*(-2d0*ERI_AO(cbas,fbas,abas,ebas)+ERI_AO(fbas,cbas,abas,ebas))
      enddo
     enddo
    enddo
   enddo
  enddo

  Ainter=czero; Binter=czero; Cinter=czero;

  ! Sigma_c_ab = \sum_cdefgh Gcd Geg Ghf v_hbdg ( v_cfae - v_fcae )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(bbas,hbas,cbas,gbas)=Ainter(bbas,hbas,cbas,gbas)+G_ao1(cbas,dbas)*ERI_AO(hbas,bbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(bbas,hbas,cbas,ebas)=Binter(bbas,hbas,cbas,ebas)+G_ao2(ebas,gbas)*Ainter(bbas,hbas,cbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do bbas=1,nBas
    do hbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(bbas,fbas,cbas,ebas)=Cinter(bbas,fbas,cbas,ebas)+G_ao3(hbas,fbas)*Binter(bbas,hbas,cbas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(bbas,fbas,cbas,ebas)*(ERI_AO(cbas,fbas,abas,ebas)-ERI_AO(fbas,cbas,abas,ebas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine 

subroutine Sigma_c_GF2B_hh_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c)

! Restricted scGF2B hh prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfh Gge v_hbgd (-2 v_aecf + v_aefc )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(hbas,bbas,gbas,cbas)=Ainter(hbas,bbas,gbas,cbas)+G_ao1(cbas,dbas)*ERI_AO(hbas,bbas,gbas,dbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(fbas,bbas,gbas,cbas)=Binter(fbas,bbas,gbas,cbas)+G_ao2(fbas,hbas)*Ainter(hbas,bbas,gbas,cbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do bbas=1,nBas
    do gbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(fbas,bbas,ebas,cbas)=Cinter(fbas,bbas,ebas,cbas)+G_ao3(gbas,ebas)*Binter(fbas,bbas,gbas,cbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(fbas,bbas,ebas,cbas)*(-2d0*ERI_AO(abas,ebas,cbas,fbas)+ERI_AO(abas,ebas,fbas,cbas))
      enddo
     enddo
    enddo
   enddo
  enddo

  Ainter=czero; Binter=czero; Cinter=czero;

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfh Gge v_hbdg ( v_aecf - v_aefc )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(hbas,bbas,gbas,cbas)=Ainter(hbas,bbas,gbas,cbas)+G_ao1(cbas,dbas)*ERI_AO(hbas,bbas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(fbas,bbas,gbas,cbas)=Binter(fbas,bbas,gbas,cbas)+G_ao2(fbas,hbas)*Ainter(hbas,bbas,gbas,cbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do fbas=1,nBas
   do bbas=1,nBas
    do gbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(fbas,bbas,ebas,cbas)=Cinter(fbas,bbas,ebas,cbas)+G_ao3(gbas,ebas)*Binter(fbas,bbas,gbas,cbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(fbas,bbas,ebas,cbas)*(ERI_AO(abas,ebas,cbas,fbas)-ERI_AO(abas,ebas,fbas,cbas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine 

subroutine Sigma_c_GF2B_ee_prime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c)

! Restricted scGF2B ee prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Geg Ghf v_gdhb (-2 v_cfae + v_fcae )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(gbas,cbas,hbas,bbas)=Ainter(gbas,cbas,hbas,bbas)+G_ao1(cbas,dbas)*ERI_AO(gbas,dbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(ebas,cbas,hbas,bbas)=Binter(ebas,cbas,hbas,bbas)+G_ao2(ebas,gbas)*Ainter(gbas,cbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do hbas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(ebas,cbas,fbas,bbas)=Cinter(ebas,cbas,fbas,bbas)+G_ao3(hbas,fbas)*Binter(ebas,cbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(ebas,cbas,fbas,bbas)*(-2d0*ERI_AO(cbas,fbas,abas,ebas)+ERI_AO(fbas,cbas,abas,ebas))
      enddo
     enddo
    enddo
   enddo
  enddo

  Ainter=czero; Binter=czero; Cinter=czero;

  ! Sigma_c_ab = \sum_cdefgh Gcd Geg Ghf v_dghb ( v_cfae - v_fcae )
  do cbas=1,nBas
   do gbas=1,nBas
    do bbas=1,nBas
     do hbas=1,nBas
      do dbas=1,nBas
       Ainter(gbas,cbas,hbas,bbas)=Ainter(gbas,cbas,hbas,bbas)+G_ao1(cbas,dbas)*ERI_AO(dbas,gbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do ebas=1,nBas
   do cbas=1,nBas
    do gbas=1,nBas
     do bbas=1,nBas
      do hbas=1,nBas
       Binter(ebas,cbas,hbas,bbas)=Binter(ebas,cbas,hbas,bbas)+G_ao2(ebas,gbas)*Ainter(gbas,cbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do hbas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Cinter(ebas,cbas,fbas,bbas)=Cinter(ebas,cbas,fbas,bbas)+G_ao3(hbas,fbas)*Binter(ebas,cbas,hbas,bbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do fbas=1,nBas
     do cbas=1,nBas
      do ebas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(ebas,cbas,fbas,bbas)*(ERI_AO(cbas,fbas,abas,ebas)-ERI_AO(fbas,cbas,abas,ebas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine 

subroutine Sigma_c_GF2B_hh_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c)

! Restricted scGF2B hh prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Gfg Ghe v_aecf (2 v_bhgd - v_hbdg )
  do abas=1,nBas
   do ebas=1,nBas
    do cbas=1,nBas
     do fbas=1,nBas
      do dbas=1,nBas
       Ainter(abas,ebas,dbas,fbas)=Ainter(abas,ebas,dbas,fbas)+G_ao1(cbas,dbas)*ERI_AO(abas,ebas,cbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do ebas=1,nBas
    do dbas=1,nBas
     do fbas=1,nBas
      do gbas=1,nBas
       Binter(abas,ebas,dbas,gbas)=Binter(abas,ebas,dbas,gbas)+G_ao2(fbas,gbas)*Ainter(abas,ebas,dbas,fbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do ebas=1,nBas
    do dbas=1,nBas
     do gbas=1,nBas
      do hbas=1,nBas
       Cinter(abas,hbas,dbas,gbas)=Cinter(abas,hbas,dbas,gbas)+G_ao3(hbas,ebas)*Binter(abas,ebas,dbas,gbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do hbas=1,nBas
     do dbas=1,nBas
      do gbas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(abas,hbas,dbas,gbas)*(2d0*ERI_AO(bbas,hbas,gbas,dbas)-ERI_AO(hbas,bbas,dbas,gbas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine

subroutine Sigma_c_GF2B_ee_primeprime(nBas,Ainter,Binter,Cinter,G_ao1,G_ao2,G_ao3,ERI_AO,Sigma_c)

! Restricted scGF2B ee prime prime

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_ao1(nBas,nBas)
  complex*16,intent(in)         :: G_ao2(nBas,nBas)
  complex*16,intent(in)         :: G_ao3(nBas,nBas)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

! Output variables
  complex*16,intent(inout)      :: Sigma_c(nBas,nBas)
  complex*16,intent(inout)      :: Ainter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Binter(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: Cinter(nBas,nBas,nBas,nBas)

  ! Sigma_c_ab = \sum_cdefgh Gcd Geh Ggf v_cfae (2 v_dgbh - v_dghb )
  do cbas=1,nBas
   do fbas=1,nBas
    do abas=1,nBas
     do ebas=1,nBas
      do dbas=1,nBas
       Ainter(dbas,fbas,abas,ebas)=Ainter(dbas,fbas,abas,ebas)+G_ao1(cbas,dbas)*ERI_AO(cbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas
   do fbas=1,nBas
    do abas=1,nBas
     do ebas=1,nBas
      do hbas=1,nBas
       Binter(dbas,fbas,abas,hbas)=Binter(dbas,fbas,abas,hbas)+G_ao2(ebas,hbas)*Ainter(dbas,fbas,abas,ebas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do dbas=1,nBas
   do fbas=1,nBas
    do abas=1,nBas
     do hbas=1,nBas
      do gbas=1,nBas
       Cinter(dbas,gbas,abas,hbas)=Cinter(dbas,gbas,abas,hbas)+G_ao3(gbas,fbas)*Binter(dbas,fbas,abas,hbas)
      enddo
     enddo
    enddo
   enddo
  enddo
  do abas=1,nBas
   do bbas=1,nBas
    do dbas=1,nBas
     do gbas=1,nBas
      do hbas=1,nBas
       Sigma_c(abas,bbas) =Sigma_c(abas,bbas)+Cinter(dbas,gbas,abas,hbas)*(2d0*ERI_AO(dbas,gbas,bbas,hbas)-ERI_AO(dbas,gbas,hbas,bbas))
      enddo
     enddo
    enddo
   enddo
  enddo

end subroutine

subroutine Sigma_c_GF2B_brut(nBas,nBas_twice,G_plus,G_minus,ERI_AO,Sigma_c_plus,Sigma_c_minus) 

! Restricted scGF2B all block M^8

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas_twice

  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

  complex*16,intent(in)         :: G_plus(nBas_twice,nBas_twice)
  complex*16,intent(in)         :: G_minus(nBas_twice,nBas_twice)

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
  complex*16,intent(inout)      :: Sigma_c_plus(nBas_twice,nBas_twice)
  complex*16,intent(inout)      :: Sigma_c_minus(nBas_twice,nBas_twice)

  ! Allocate arrays

  allocate(G_he_plus(nBas,nBas))
  allocate(G_hh_plus(nBas,nBas))
  allocate(G_ee_plus(nBas,nBas))
  allocate(G_eh_plus(nBas,nBas))
  allocate(G_he_minus(nBas,nBas))
  allocate(G_hh_minus(nBas,nBas))
  allocate(G_ee_minus(nBas,nBas))
  allocate(G_eh_minus(nBas,nBas))
  allocate(Sigma_he_plus(nBas,nBas)) 
  allocate(Sigma_hh_plus(nBas,nBas))
  allocate(Sigma_ee_plus(nBas,nBas))
  allocate(Sigma_eh_plus(nBas,nBas))
  allocate(Sigma_he_minus(nBas,nBas))
  allocate(Sigma_hh_minus(nBas,nBas))
  allocate(Sigma_ee_minus(nBas,nBas))
  allocate(Sigma_eh_minus(nBas,nBas))

  ! Initialize 
  
  G_he_plus(1:nBas,1:nBas)=G_plus(           1:nBas,           1:nBas)
  G_hh_plus(1:nBas,1:nBas)=G_plus(           1:nBas,nBas+1:nBas_twice)    
  G_ee_plus(1:nBas,1:nBas)=G_plus(nBas+1:nBas_twice,           1:nBas)   
  G_eh_plus(1:nBas,1:nBas)=G_plus(nBas+1:nBas_twice,nBas+1:nBas_twice)  
  G_he_minus(1:nBas,1:nBas)=G_minus(           1:nBas,           1:nBas)
  G_hh_minus(1:nBas,1:nBas)=G_minus(           1:nBas,nBas+1:nBas_twice)    
  G_ee_minus(1:nBas,1:nBas)=G_minus(nBas+1:nBas_twice,           1:nBas)   
  G_eh_minus(1:nBas,1:nBas)=G_minus(nBas+1:nBas_twice,nBas+1:nBas_twice)  
  Sigma_he_plus=czero
  Sigma_hh_plus=czero
  Sigma_ee_plus=czero
  Sigma_eh_plus=czero
  Sigma_he_minus=czero
  Sigma_hh_minus=czero
  Sigma_ee_minus=czero
  Sigma_eh_minus=czero

  ! M^8 loop
  write(*,*)
  write(*,*) 'Computing Bog. Sigma_c (M^8)'
  write(*,*)
  do abas=1,nBas
   do bbas=1,nBas
    do cbas=1,nBas
     do dbas=1,nBas
      do ebas=1,nBas
       do fbas=1,nBas
        do gbas=1,nBas
         do hbas=1,nBas
           ! he 2'
           Integ_val=ERI_AO(abas,ebas,cbas,fbas)*(2d0*ERI_AO(dbas,gbas,bbas,hbas)-ERI_AO(dbas,gbas,hbas,bbas))
           Sigma_he_plus(abas,bbas) = Sigma_he_plus(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_he_plus(fbas,gbas) *G_he_minus(hbas,ebas)            
           Sigma_he_minus(abas,bbas)=Sigma_he_minus(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_he_minus(fbas,gbas)*G_he_plus(hbas,ebas)            
           ! he 2''
           Integ_val=ERI_AO(abas,ebas,cbas,fbas)*(ERI_AO(dbas,gbas,hbas,bbas)-2d0*ERI_AO(dbas,gbas,bbas,hbas)) &
                    +ERI_AO(abas,ebas,fbas,cbas)*(ERI_AO(dbas,gbas,bbas,hbas)-ERI_AO(dbas,gbas,hbas,bbas))
           Sigma_he_plus(abas,bbas) = Sigma_he_plus(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_hh_plus(fbas,hbas) *G_ee_minus(gbas,ebas)            
           Sigma_he_minus(abas,bbas)=Sigma_he_minus(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_hh_minus(fbas,hbas)*G_ee_plus(gbas,ebas) 
           ! eh 2'
           Integ_val=ERI_AO(cbas,fbas,abas,ebas)*(2d0*ERI_AO(bbas,hbas,dbas,gbas)-ERI_AO(bbas,hbas,gbas,dbas))
           Sigma_eh_plus(abas,bbas) = Sigma_eh_plus(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_eh_plus(fbas,gbas) *G_eh_minus(hbas,ebas)            
           Sigma_eh_minus(abas,bbas)=Sigma_eh_minus(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_eh_minus(fbas,gbas)*G_eh_plus(hbas,ebas)            
           ! eh 2''
           Integ_val=ERI_AO(cbas,fbas,abas,ebas)*(ERI_AO(bbas,hbas,gbas,dbas)-2d0*ERI_AO(bbas,hbas,dbas,gbas)) &
                    +ERI_AO(cbas,fbas,ebas,abas)*(ERI_AO(bbas,hbas,dbas,gbas)-ERI_AO(bbas,hbas,gbas,dbas))
           Sigma_eh_plus(abas,bbas) = Sigma_eh_plus(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_hh_plus(ebas,gbas) *G_ee_minus(hbas,fbas)            
           Sigma_eh_minus(abas,bbas)=Sigma_eh_minus(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_hh_minus(ebas,gbas)*G_ee_plus(hbas,fbas) 
           ! hh 2'
           Integ_val=ERI_AO(abas,ebas,cbas,fbas)*(ERI_AO(hbas,bbas,dbas,gbas)-2d0*ERI_AO(hbas,bbas,gbas,dbas)) &
                    +ERI_AO(abas,ebas,fbas,cbas)*(ERI_AO(hbas,bbas,gbas,dbas)-ERI_AO(hbas,bbas,dbas,gbas))
           Sigma_hh_plus(abas,bbas) = Sigma_hh_plus(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_he_plus(fbas,hbas) *G_he_minus(gbas,ebas)            
           Sigma_hh_minus(abas,bbas)=Sigma_hh_minus(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_he_minus(fbas,hbas)*G_he_plus(gbas,ebas) 
           ! hh 2''
           Integ_val=ERI_AO(abas,ebas,cbas,fbas)*(2d0*ERI_AO(hbas,bbas,gbas,dbas)-ERI_AO(hbas,bbas,dbas,gbas))
           Sigma_hh_plus(abas,bbas) = Sigma_hh_plus(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_hh_plus(fbas,gbas) *G_ee_minus(hbas,ebas)            
           Sigma_hh_minus(abas,bbas)=Sigma_hh_minus(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_hh_minus(fbas,gbas)*G_ee_plus(hbas,ebas)            
           ! ee 2'
           Integ_val=ERI_AO(cbas,fbas,abas,ebas)*(ERI_AO(gbas,dbas,bbas,hbas)-2d0*ERI_AO(gbas,dbas,hbas,bbas)) &
                    +ERI_AO(cbas,fbas,ebas,abas)*(ERI_AO(gbas,dbas,hbas,bbas)-ERI_AO(gbas,dbas,bbas,hbas))
           Sigma_ee_plus(abas,bbas) = Sigma_ee_plus(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_he_plus(ebas,gbas) *G_he_minus(hbas,fbas)            
           Sigma_ee_minus(abas,bbas)=Sigma_ee_minus(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_he_minus(ebas,gbas)*G_he_plus(hbas,fbas) 
           ! ee 2''
           Integ_val=ERI_AO(cbas,fbas,abas,ebas)*(2d0*ERI_AO(gbas,dbas,hbas,bbas)-ERI_AO(gbas,dbas,bbas,hbas))
           Sigma_ee_plus(abas,bbas) = Sigma_ee_plus(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_hh_plus(ebas,hbas) *G_ee_minus(gbas,fbas)            
           Sigma_ee_minus(abas,bbas)=Sigma_ee_minus(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_hh_minus(ebas,hbas)*G_ee_plus(gbas,fbas)            
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo

  ! Set Sigma_c Gorkov

  Sigma_c_plus(           1:nBas,           1:nBas)=Sigma_he_plus(1:nBas,1:nBas)
  Sigma_c_plus(           1:nBas,nBas+1:nBas_twice)=Sigma_hh_plus(1:nBas,1:nBas)    
  Sigma_c_plus(nBas+1:nBas_twice,           1:nBas)=Sigma_ee_plus(1:nBas,1:nBas)   
  Sigma_c_plus(nBas+1:nBas_twice,nBas+1:nBas_twice)=Sigma_eh_plus(1:nBas,1:nBas)  
  Sigma_c_minus(           1:nBas,           1:nBas)=Sigma_he_minus(1:nBas,1:nBas)
  Sigma_c_minus(           1:nBas,nBas+1:nBas_twice)=Sigma_hh_minus(1:nBas,1:nBas)    
  Sigma_c_minus(nBas+1:nBas_twice,           1:nBas)=Sigma_ee_minus(1:nBas,1:nBas)   
  Sigma_c_minus(nBas+1:nBas_twice,nBas+1:nBas_twice)=Sigma_eh_minus(1:nBas,1:nBas)  

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

