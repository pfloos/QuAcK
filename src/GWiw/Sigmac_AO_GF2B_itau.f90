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



