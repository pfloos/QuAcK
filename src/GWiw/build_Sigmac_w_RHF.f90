subroutine build_Sigmac_w_RHF(nOrb,nO,nE,verbose,wtest,eHF,nfreqs,ntimes,wweight,wcoord,vMAT,&
                              Sigma_c_mo)

! Restricted Sigma_c(E)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nE
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: vMAT(nOrb*nOrb,nOrb*nOrb)

  complex *16,intent(in)        :: wtest(nE)

! Local variables

  integer                       :: ifreq,iE
  integer                       :: iorb,jorb,korb,lorb,porb
  integer                       :: nOrb2

  double precision              :: chem_pot
  double precision              :: eta
  double precision,external     :: Heaviside_step
  double precision,allocatable  :: Tmp_mo_w(:,:)

  complex *16                   :: weval
  complex *16,allocatable       :: Chi0_mo_w(:,:)
  complex *16,allocatable       :: G_mo_1(:,:)
  complex *16,allocatable       :: G_mo_2(:,:)

! Ouput variables

  complex *16,intent(out)       :: Sigma_c_mo(nE,nOrb,nOrb)

!

  nOrb2=nOrb*nOrb
  Sigma_c_mo=czero
  allocate(Chi0_mo_w(nOrb2,nOrb2),Tmp_mo_w(nOrb2,nOrb2))
  allocate(G_mo_1(nOrb,nOrb),G_mo_2(nOrb,nOrb))

! Imaginary freqs contribution

  do ifreq=1,nfreqs
   
   ! Initialize
   Tmp_mo_w=0d0
   do iorb=1,nOrb2
    Tmp_mo_w(iorb,iorb)=1d0  
   enddo

   ! Xo (iw)
   if(ntimes>0) then
    call Gitau2Chi0iw_mo_RHF(nOrb,nO,0,eHF,ntimes,wcoord(ifreq),Chi0_mo_w)
   else
    eta=0d0
    call Xoiw_mo_RHF(nOrb,nO,eta,eHF,im*wcoord(ifreq),Chi0_mo_w)
   endif

   ! Xo (iw) -> Wp (iw)
   Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
   call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
   Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
   Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
   Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp(iw) in MO

   ! Use Wp (iw) to build all Sigma_c(E)
    eta=0d0
    do iE=1,nE
     ! Build G(iw+wtest)
     weval=wtest(iE)+im*wcoord(ifreq)
     call G_MO_RHF(nOrb,nO,eta,eHF,weval,G_mo_1)
     weval=wtest(iE)-im*wcoord(ifreq)
     call G_MO_RHF(nOrb,nO,eta,eHF,weval,G_mo_2)
     ! Sigma_c(E)
     do iorb=1,nOrb
      do jorb=1,nOrb
       do korb=1,nOrb
        Sigma_c_mo(iE,iorb,jorb)=Sigma_c_mo(iE,iorb,jorb)-(G_mo_1(korb,korb)+G_mo_2(korb,korb))  &
                             *Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(korb-1)*nOrb)        &
                             *wweight(ifreq)/(2d0*pi)
       enddo
      enddo
     enddo
    enddo

  enddo

! Residues contributions

  do iE=1,nE

   ! Contour deformation residues
   if(abs(aimag(wtest(iE)))<1e-12) then ! wtest is real and we may have to add residues contributions

     eta=0.00001

     ! Align the poles of G
     chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
     eHF(:) = eHF(:)-chem_pot

     ! Occupied residues
     do porb=1,nO
      if(Heaviside_step(eHF(porb)-Real(wtest(iE)))>0d0) then
       ! Building Wp(E)
       Tmp_mo_w=0d0
       do iorb=1,nOrb2
        Tmp_mo_w(iorb,iorb)=1d0
       enddo
       weval=eHF(porb)-Real(wtest(iE))
       call Xoiw_mo_RHF(nOrb,nO,eta,eHF,weval,Chi0_mo_w)
       Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
       call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
       Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp in MO
       do iorb=1,nOrb
        do jorb=1,nOrb
         do korb=1,nO
          Sigma_c_mo(iE,iorb,jorb)=Sigma_c_mo(iE,iorb,jorb) &
                               -Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(korb-1)*nOrb)
         enddo
        enddo
       enddo
      endif
     enddo

     ! Virtual residues
     do porb=nO+1,nOrb
      if(Heaviside_step(Real(wtest(iE))-eHF(porb))>0d0) then
       ! Building Wp(E)
       Tmp_mo_w=0d0
       do iorb=1,nOrb2
        Tmp_mo_w(iorb,iorb)=1d0
       enddo
       weval=Real(wtest(iE))-eHF(porb)
       call Xoiw_mo_RHF(nOrb,nO,eta,eHF,weval,Chi0_mo_w)
       Tmp_mo_w(:,:)=Tmp_mo_w(:,:)-matmul(Real(Chi0_mo_w(:,:)),vMAT(:,:))
       call inverse_matrix(nOrb2,Tmp_mo_w,Tmp_mo_w)
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),Real(Chi0_mo_w(:,:)))
       Tmp_mo_w(:,:)=matmul(Tmp_mo_w(:,:),vMAT(:,:))
       Tmp_mo_w(:,:)=matmul(vMAT(:,:),Tmp_mo_w(:,:)) ! Now Tmp_mo_w is Wp in MO
       do iorb=1,nOrb
        do jorb=1,nOrb
         do korb=nO+1,nOrb
          Sigma_c_mo(iE,iorb,jorb)=Sigma_c_mo(iE,iorb,jorb) &
                               +Tmp_mo_w(1+(korb-1)+(iorb-1)*nOrb,1+(jorb-1)+(korb-1)*nOrb)
         enddo
        enddo
       enddo
      endif
     enddo

     ! Recover eHF
     eHF(:) = eHF(:)+chem_pot

   endif

  enddo


! Print results
 
  if(verbose/=0) then 
   do iE=1,nE
    write(*,*) ' ' 
    write(*,'(a,f15.8,a,f15.8,a)') ' RHF Sigma_c(wtest) in MO for wtest=(',Real(wtest(iE)), &
    ",",Aimag(wtest(iE)),")"
    write(*,*) ' ' 
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Real(Sigma_c_mo(iE,iorb,:))
    enddo
    write(*,*) ' ' 
   enddo
  endif

  ! Deallocate arrays
  deallocate(Chi0_mo_w)
  deallocate(Tmp_mo_w)
  deallocate(G_mo_1,G_mo_2)

end subroutine

