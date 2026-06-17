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

