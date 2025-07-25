subroutine build_Sigmac_w_RHFB(nOrb,nOrb_twice,nE,eta,verbose,wtest,eHFB,nfreqs,ntimes,&
                              wweight,wcoord,vMAT,U_QP,Sigma_he_c_mo,Sigma_hh_c_mo,    &
                              Sigma_eh_c_mo,Sigma_ee_c_mo,doSigc_eh,doSigc_ee)

! Restricted Sigma_c(E)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doSigc_eh
  logical,intent(in)            :: doSigc_ee

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nE
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: eta
  double precision,intent(inout):: eHFB(nOrb_twice)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)

  complex *16,intent(in)        :: wtest(nE)

! Local variables

  integer                       :: ifreq,iE
  integer                       :: iorb,jorb,korb,lorb,Istate
  integer                       :: nOrb2

  double precision              :: chem_pot
  double precision,external     :: Heaviside_step
  double precision,allocatable  :: Tmp_mo_w(:,:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)

  complex *16                   :: weval
  complex *16,allocatable       :: Chi0_mo_w(:,:)
  complex *16,allocatable       :: G_mo_1(:,:)
  complex *16,allocatable       :: G_mo_2(:,:)
  complex *16,allocatable       :: G_mo_3(:,:)
  complex *16,allocatable       :: G_mo_4(:,:)
  complex *16,allocatable       :: G_mo_5(:,:)
  complex *16,allocatable       :: G_mo_6(:,:)
  complex *16,allocatable       :: G_mo_7(:,:)
  complex *16,allocatable       :: G_mo_8(:,:)

! Ouput variables

  complex *16,intent(out)       :: Sigma_he_c_mo(nE,nOrb,nOrb)
  complex *16,intent(out)       :: Sigma_hh_c_mo(nE,nOrb,nOrb)
  complex *16,intent(out)       :: Sigma_eh_c_mo(nE,nOrb,nOrb)
  complex *16,intent(out)       :: Sigma_ee_c_mo(nE,nOrb,nOrb)

!

  nOrb2=nOrb*nOrb
  Sigma_he_c_mo=czero
  Sigma_hh_c_mo=czero
  Sigma_eh_c_mo=czero
  Sigma_ee_c_mo=czero
  allocate(Chi0_mo_w(nOrb2,nOrb2),Tmp_mo_w(nOrb2,nOrb2))
  allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb))
  allocate(G_mo_1(nOrb,nOrb),G_mo_2(nOrb,nOrb))
  allocate(G_mo_3(nOrb,nOrb),G_mo_4(nOrb,nOrb))
  allocate(G_mo_5(nOrb,nOrb),G_mo_6(nOrb,nOrb))
  allocate(G_mo_7(nOrb,nOrb),G_mo_8(nOrb,nOrb))
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

! Imaginary freqs contribution

  do ifreq=1,nfreqs
   
   ! Initialize
   Tmp_mo_w=0d0
   do iorb=1,nOrb2
    Tmp_mo_w(iorb,iorb)=1d0  
   enddo

   ! Xo (iw)
   if(ntimes>0) then
    call Gitau2Chi0iw_mo_RHFB(nOrb,nOrb_twice,0,eHFB,ntimes,wcoord(ifreq),U_QP,Chi0_mo_w)
   else
    call Xoiw_mo_RHFB(nOrb,nOrb_twice,eta,eHFB,im*wcoord(ifreq),Mat1,Mat2,Chi0_mo_w)
   endif

   ! Xo (iw) -> Wp (iw)
   Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
   call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
   Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
   Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
   Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp(iw) in MO

   ! Use Wp (iw) to build all Sigma_c(E)
    do iE=1,nE
     ! Build G(iw+wtest)
     ! Ghe
      weval=wtest(iE)+im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat1,Mat2,Mat2,G_mo_1)
      weval=wtest(iE)-im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat1,Mat2,Mat2,G_mo_2)
     ! Ghh
      weval=wtest(iE)+im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,-Mat2,Mat1,G_mo_3)
      weval=wtest(iE)-im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,-Mat2,Mat1,G_mo_4)
     ! Geh
     if(doSigc_eh) then
      weval=wtest(iE)+im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat2,Mat2,Mat1,Mat1,G_mo_5)
      weval=wtest(iE)-im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat2,Mat2,Mat1,Mat1,G_mo_6)
     endif
     ! Gee
     if(doSigc_ee) then
      weval=wtest(iE)+im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat2,Mat1,Mat1,-Mat2,G_mo_7)
      weval=wtest(iE)-im*wcoord(ifreq)
      call G_MO_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat2,Mat1,Mat1,-Mat2,G_mo_8)
     endif
     ! Sigma_c(E)
     do iorb=1,nOrb
      do jorb=1,nOrb
       do korb=1,nOrb
        do lorb=1,nOrb
         ! This is  G^he W
         Sigma_he_c_mo(iE,iorb,jorb)=Sigma_he_c_mo(iE,iorb,jorb)-(G_mo_1(korb,lorb)+G_mo_2(korb,lorb))  &
                                    *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)        &
                                    *wweight(ifreq)/(2d0*pi)
         ! This is -G^hh W
         Sigma_hh_c_mo(iE,iorb,jorb)=Sigma_hh_c_mo(iE,iorb,jorb)+(G_mo_3(korb,lorb)+G_mo_4(korb,lorb))  &
                                    *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)        &
                                    *wweight(ifreq)/(2d0*pi)
         ! This is  G^eh W
         if(doSigc_eh) then
         Sigma_eh_c_mo(iE,iorb,jorb)=Sigma_eh_c_mo(iE,iorb,jorb)-(G_mo_5(korb,lorb)+G_mo_6(korb,lorb))  &
                                    *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)        &
                                    *wweight(ifreq)/(2d0*pi)
         endif
         if(doSigc_ee) then
         ! This is -G^ee W
         Sigma_ee_c_mo(iE,iorb,jorb)=Sigma_ee_c_mo(iE,iorb,jorb)+(G_mo_7(korb,lorb)+G_mo_8(korb,lorb))  &
                                    *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)        &
                                    *wweight(ifreq)/(2d0*pi)
         endif
        enddo
       enddo
      enddo
     enddo
    enddo

  enddo

! Residues contributions

  do iE=1,nE

   ! Contour deformation residues
   if(abs(aimag(wtest(iE)))<1e-12) then ! wtest is real and we may have to add residues contributions

     ! In general, we only need NEGATIVE REAL (OR PURELY IMAGINARY) wtest values when building  Sigma_c^Gorkov
     ! because for POSITIVE REAL wtest values, we have the following equalities
     !  Sigma_c^he(wtest) = -Sigma_c^eh(-wtest) 
     !  Sigma_c^eh(wtest) = -Sigma_c^he(-wtest) 
     !  Sigma_c^hh(wtest) = Sigma_c^ee(-wtest) 
     !  Sigma_c^ee(wtest) = Sigma_c^hh(-wtest) 
     do Istate=1,nOrb
      if(Heaviside_step(eHFB(Istate)-Real(wtest(iE)))>0d0) then
       ! Building Wp(E)
       Tmp_mo_w=0d0
       do iorb=1,nOrb2
        Tmp_mo_w(iorb,iorb)=1d0
       enddo
       weval=eHFB(Istate)-Real(wtest(iE))
       call Xoiw_mo_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,Chi0_mo_w)
       Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
       call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
       Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp in MO
       do iorb=1,nOrb
        do jorb=1,nOrb
         do korb=1,nOrb
          do lorb=1,nOrb
           !  G^he W
           Sigma_he_c_mo(iE,iorb,jorb)=Sigma_he_c_mo(iE,iorb,jorb)-Mat1(korb,Istate)*Mat1(lorb,Istate) &
                                      *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)
           ! -G^hh W
           Sigma_hh_c_mo(iE,iorb,jorb)=Sigma_hh_c_mo(iE,iorb,jorb)+Mat1(korb,Istate)*Mat2(lorb,Istate) &
                                      *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)
           !  G^eh W
           if(doSigc_eh) then
           Sigma_eh_c_mo(iE,iorb,jorb)=Sigma_eh_c_mo(iE,iorb,jorb)-Mat2(korb,Istate)*Mat2(lorb,Istate) &
                                      *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)
           endif
           ! -G^ee W
           if(doSigc_ee) then
           Sigma_ee_c_mo(iE,iorb,jorb)=Sigma_ee_c_mo(iE,iorb,jorb)+Mat2(korb,Istate)*Mat1(lorb,Istate) &
                                      *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)
           endif
          enddo
         enddo
        enddo
       enddo
      endif
      if(Heaviside_step(Real(wtest(iE))-abs(eHFB(Istate)))>0d0) then
       ! Building Wp(E)
       Tmp_mo_w=0d0
       do iorb=1,nOrb2
        Tmp_mo_w(iorb,iorb)=1d0
       enddo
       weval=-(Real(wtest(iE))-abs(eHFB(Istate)))
       call Xoiw_mo_RHFB(nOrb,nOrb_twice,eta,eHFB,weval,Mat1,Mat2,Chi0_mo_w)
       Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
       call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
       Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp in MO
       do iorb=1,nOrb
        do jorb=1,nOrb
         do korb=1,nOrb
          do lorb=1,nOrb
           !  G^he W
           Sigma_he_c_mo(iE,iorb,jorb)=Sigma_he_c_mo(iE,iorb,jorb)+Mat2(korb,Istate)*Mat2(lorb,Istate) &
                                      *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)
           ! -G^hh W
           Sigma_hh_c_mo(iE,iorb,jorb)=Sigma_hh_c_mo(iE,iorb,jorb)+Mat2(korb,Istate)*Mat1(lorb,Istate) &
                                      *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)
           !  G^eh W
           if(doSigc_eh) then
           Sigma_eh_c_mo(iE,iorb,jorb)=Sigma_eh_c_mo(iE,iorb,jorb)+Mat1(korb,Istate)*Mat1(lorb,Istate) &
                                      *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)
           endif
           ! -G^ee W
           if(doSigc_ee) then
           Sigma_ee_c_mo(iE,iorb,jorb)=Sigma_ee_c_mo(iE,iorb,jorb)+Mat1(korb,Istate)*Mat2(lorb,Istate) &
                                      *Tmp_mo_w(1+(iorb-1)+(korb-1)*nOrb,1+(jorb-1)+(lorb-1)*nOrb)
           endif
          enddo
         enddo
        enddo
       enddo
      endif
     enddo

   endif

  enddo

! Print results

  if(verbose/=0) then  
   do iE=1,nE
    write(*,*) ' ' 
    write(*,'(a,f15.8,a,f15.8,a)') ' HFB Sigma_he/hh_c(wtest) in MO for wtest=(',Real(wtest(iE)), &
    ",",Aimag(wtest(iE)),")"
    write(*,*) ' ' 
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Real(Sigma_he_c_mo(iE,iorb,:))
    enddo
    write(*,*) ' ' 
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Real(Sigma_hh_c_mo(iE,iorb,:))
    enddo
    write(*,*) ' ' 
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Real(Sigma_eh_c_mo(iE,iorb,:))
    enddo
    write(*,*) ' ' 
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Real(Sigma_ee_c_mo(iE,iorb,:))
    enddo
   enddo
  endif

  ! Deallocate arrays
  deallocate(Chi0_mo_w)
  deallocate(Tmp_mo_w)
  deallocate(G_mo_1,G_mo_2)
  deallocate(G_mo_3,G_mo_4)
  deallocate(G_mo_5,G_mo_6)
  deallocate(G_mo_7,G_mo_8)
  deallocate(Mat1,Mat2)

end subroutine

