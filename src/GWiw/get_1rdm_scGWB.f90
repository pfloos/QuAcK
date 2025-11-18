subroutine get_1rdm_scGWB(nBas,nBas_twice,nfreqs,chem_pot,S,H_ao_hfb,Sigma_c_w_ao,wcoord,wweight, &
                          G_ao,G_ao_iw_hfb,DeltaG_ao_iw,R_ao,R_ao_hfb,trace_1_rdm) 

! Compute the scGW 1RDM

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas_twice
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: chem_pot
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: H_ao_hfb(nBas_twice,nBas_twice)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: R_ao_hfb(nBas_twice,nBas_twice)
  complex*16,intent(in)         :: G_ao_iw_hfb(nfreqs,nBas_twice,nBas_twice)
  complex*16,intent(in)         :: Sigma_c_w_ao(nfreqs,nBas_twice,nBas_twice)

! Local variables

  integer                       :: ifreq
  integer                       :: ibas,jbas
  double precision,allocatable  :: Stmp(:,:)
  double precision,allocatable  :: Stmp2(:,:)
  complex*16                    :: weval_cpx

! Output variables

  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: R_ao(nBas_twice,nBas_twice)
  complex*16,intent(out)        :: G_ao(nBas_twice,nBas_twice)
  complex*16,intent(out)        :: DeltaG_ao_iw(nfreqs,nBas_twice,nBas_twice)

  allocate(Stmp(nBas_twice,nBas_twice))
  allocate(Stmp2(nBas_twice,nBas_twice))
  Stmp=0d0
  Stmp2=0d0
  Stmp(1:nBas           ,1:nBas           ) =  S(1:nBas,1:nBas)
  Stmp(nBas+1:nBas_twice,nBas+1:nBas_twice) = -S(1:nBas,1:nBas) ! The minus will change the sign of the chem. pot.
  Stmp2(1:nBas           ,1:nBas           ) = S(1:nBas,1:nBas)
  Stmp2(nBas+1:nBas_twice,nBas+1:nBas_twice) = S(1:nBas,1:nBas)
  R_ao=0d0
  DeltaG_ao_iw=czero
  do ifreq=1,nfreqs
   weval_cpx=im*wcoord(ifreq)
   ! Setting G(w) = [ (w+chem_pot)S - F - Sigma_c(w) ]^-1
   G_ao(:,:)= weval_cpx*Stmp2(:,:) + chem_pot*Stmp(:,:) - H_ao_hfb(:,:) - Sigma_c_w_ao(ifreq,:,:) ! G(iw)^-1
   call complex_inverse_matrix(nBas_twice,G_ao,G_ao)                                              ! G(iw)
   G_ao(:,:)=G_ao(:,:)-G_ao_iw_hfb(ifreq,:,:)                                                     ! G_corr(iw) = G(iw) - Go(iw) 
   DeltaG_ao_iw(ifreq,:,:)=G_ao(:,:)
   R_ao(:,:) = R_ao(:,:) + wweight(ifreq)*real(G_ao(:,:))      ! P_corr = 1/(2 pi) int_-Infty ^Infty G_corr(iw) dw = 1/pi int_0 ^Infty Re[ G_corr(iw) ] dw
  enddo
  R_ao(:,:) = R_ao(:,:)/pi + R_ao_hfb(:,:)                     ! Build the correlated R_ao (only 1 spin channel)
  trace_1_rdm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+R_ao(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo
  deallocate(Stmp)
  deallocate(Stmp2)

end subroutine

