subroutine get_1rdm_scGW(nBas,nfreqs,chem_pot,S,F_ao,Sigma_c_w_ao,wcoord,wweight, &
                         G_ao,G_ao_iw_hf,DeltaG_ao_iw,P_ao,P_ao_hf,trace_1_rdm) 

! Compute the scGW 1RDM

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs
  double precision,intent(in)   :: chem_pot
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: P_ao_hf(nBas,nBas)
  complex*16,intent(in)         :: G_ao_iw_hf(nfreqs,nBas,nBas)
  complex*16,intent(in)         :: Sigma_c_w_ao(nfreqs,nBas,nBas)

! Local variables

  integer                       :: ifreq
  integer                       :: ibas,jbas
  complex*16                    :: weval_cpx

! Output variables

  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: P_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao(nBas,nBas)
  complex*16,intent(out)        :: DeltaG_ao_iw(nfreqs,nBas,nBas)

  P_ao=0d0
  DeltaG_ao_iw=czero
  do ifreq=1,nfreqs
   weval_cpx=im*wcoord(ifreq)
   ! Setting G(w) = [ (w+chem_pot)S - F - Sigma_c(w) ]^-1
   G_ao(:,:)= (weval_cpx + chem_pot)*S(:,:) - F_ao(:,:) - Sigma_c_w_ao(ifreq,:,:) ! G(iw)^-1
   call complex_inverse_matrix(nBas,G_ao,G_ao)                                    ! G(iw)
   G_ao(:,:)=G_ao(:,:)-G_ao_iw_hf(ifreq,:,:)                                      ! G_corr(iw) = G(iw) - Go(iw) 
   DeltaG_ao_iw(ifreq,:,:)=G_ao(:,:)
   P_ao(:,:) = P_ao(:,:) + wweight(ifreq)*real(G_ao(:,:))      ! P_corr = 1/(2 pi) int_-Infty ^Infty G_corr(iw) dw = 1/pi int_0 ^Infty Re[ G_corr(iw) ] dw
  enddo
  P_ao(:,:) = 2d0*P_ao(:,:)/pi + P_ao_hf(:,:)                  ! Times 2 to sum both spin channels
  trace_1_rdm=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+P_ao(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo

end subroutine

