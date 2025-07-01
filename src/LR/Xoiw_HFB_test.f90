subroutine Xoiw_HFB_tests(nBas,nOrb,nOrb_twice,cHFB,eHFB,nfreqs,ntimes,wweight,wcoord,  &
                          U_QP,ERI_AO)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: eHFB(nOrb_twice)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: fulltest

  integer                       :: ibas,jbas,lbas,kbas,ifreq
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: nBas2,nOrb2

  double precision              :: start_Xoiw   ,end_Xoiw     ,t_Xoiw
  double precision              :: EcRPA,EcGM,trace,trace2,trace3
  double precision              :: eta
  double precision,allocatable  :: Chi0_ao_iw_v(:,:)
  double precision,allocatable  :: Chi_ao_iw_v(:,:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: eigval_eps(:)
  double precision,allocatable  :: eps(:,:)
  double precision,allocatable  :: epsm1(:,:)
  double precision,allocatable  :: eigenv_eps(:,:)
  double precision,allocatable  :: vMat(:,:)
  double precision,allocatable  :: Wp_tmp(:,:)
  double precision,allocatable  :: Mat1(:,:)
  double precision,allocatable  :: Mat2(:,:)
  double precision,allocatable  :: Mat3(:,:)
  double precision,allocatable  :: Mat4(:,:)
  double precision,allocatable  :: Wp_AO(:,:,:,:)
  double precision,allocatable  :: Wp_MO(:,:,:,:)

  complex *16                   :: wtest,weval
  complex *16,allocatable       :: Sigma_he_c_ao(:,:)
  complex *16,allocatable       :: Sigma_hh_c_ao(:,:)
  complex *16,allocatable       :: G_ao_1(:,:)
  complex *16,allocatable       :: G_ao_2(:,:)
  complex *16,allocatable       :: G_ao_3(:,:)
  complex *16,allocatable       :: G_ao_4(:,:)
  complex *16,allocatable       :: cHFB_complex(:,:)
  complex *16,allocatable       :: Chi0_mo_iw(:,:)
  complex *16,allocatable       :: Chi0_ao_iw(:,:,:)
!

  fulltest=.true.      ! TODO adjust it to print Chi0(iw), Wp, and Sigma_c
  nBas2=nBas*nBas
  nOrb2=nOrb*nOrb
  wtest=0.000005967*im ! TODO use test values

!------------------------------------------------------------------------
! Build G(i tau) in AO basis
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'*******************************************'
  write(*,*)'* Use HFB Xo(i w) to build  X(i w) and    *'
  write(*,*)'*       compute EcRPA and EcGM.           *'
  write(*,*)'* Then, build Wp(i w) and Sigma_c(wtest)  *'
  write(*,*)'*******************************************'
  write(*,*)

  allocate(eps(nBas2,nBas2),epsm1(nBas2,nBas2),eigenv_eps(nBas2,nBas2),eigval_eps(nBas2))
  allocate(Chi0_ao_iw(nfreqs,nBas2,nBas2))
  allocate(Chi0_ao_iw_v(nBas2,nBas2),Chi_ao_iw_v(nBas2,nBas2))
  allocate(Wp_AO(nBas,nBas,nBas,nBas),Wp_MO(nOrb,nOrb,nOrb,nOrb),Wp_tmp(nOrb*nOrb,nOrb*nOrb))
  allocate(Wp_ao_iw(nBas2,nBas2))
  allocate(vMat(nBas2,nBas2))
  allocate(Mat1(nOrb,nOrb))
  allocate(Mat2(nOrb,nOrb))
  allocate(Mat3(nOrb,nOrb))
  allocate(Mat4(nOrb,nOrb))
  allocate(Sigma_he_c_ao(nBas,nBas),G_ao_1(nBas,nBas),G_ao_2(nBas,nBas))
  allocate(Sigma_hh_c_ao(nBas,nBas),G_ao_3(nBas,nBas),G_ao_4(nBas,nBas))
  Sigma_he_c_ao=czero
  Sigma_hh_c_ao=czero

!-----------------------------!
! Store v also as a 2D matrix !
!-----------------------------!
 
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas
      vMat(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)=ERI_AO(ibas,jbas,kbas,lbas)
     enddo
    enddo
   enddo
  enddo

!-----------------------------!
! Build X0(i w) from G(i tau) !
!-----------------------------!

  call Gitau2Chi0iw_HFB(nBas,nBas2,nOrb,nOrb_twice,cHFB,eHFB,nfreqs,ntimes,wcoord,U_QP,Chi0_ao_iw)

  if(fulltest) then

   ifreq=1; eta=0.00001; ! TODO select a frequency of wcoord
   allocate(cHFB_complex(nOrb,nOrb))
   allocate(Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb))

   Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
   Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)

   call Xoiw_HFB(nOrb,nOrb_twice,eta,eHFB,wcoord(ifreq)*im,Mat1,Mat2,Chi0_mo_iw)

   write(*,'(a,f15.8,a)') ' HFB Xo in MO for wcoord=(',wcoord(ifreq),")" 
   write(*,*) ' '
   do iorb=1,nOrb2
    write(*,'(*(f10.5))') Real(Chi0_mo_iw(iorb,:))
   enddo
   write(*,*) ' '

   deallocate(cHFB_complex)
   deallocate(Chi0_mo_iw)

  endif

!----------------------!
! Use Xo(i w) as usual !
!----------------------!

  call wall_time(start_Xoiw)

  EcRPA=0d0; EcGM=0d0;
  do ifreq=1,nfreqs

    ! Initialization
    trace=0d0; trace2=0d0; trace3=0d0;
    eigval_eps(:)=0d0
    eigenv_eps(:,:)=0d0
    eps(:,:)=0d0
   
    ! Build Xo v
    Chi0_ao_iw_v(:,:)=matmul(Real(Chi0_ao_iw(ifreq,:,:)),vMat(:,:))
   
    ! Tr [ Xo v ] and define eps
    do ibas=1,nBas2
     eps(ibas,ibas) = 1d0
     trace=trace+Chi0_ao_iw_v(ibas,ibas)
    enddo
    eps(:,:) = eps(:,:) - Chi0_ao_iw_v(:,:)
   
    ! For GW 
    ! a) build eps^-1
    ! b) Set X v = eps^-1 Xo v 
    call inverse_matrix(nBas2,eps,epsm1)
    Chi_ao_iw_v(:,:)=matmul(epsm1(:,:),Real(Chi0_ao_iw(ifreq,:,:)))
    Chi_ao_iw_v(:,:)=matmul(Chi_ao_iw_v(:,:),vMat(:,:))
   
    ! Tr [ X v ] for GM
    do ibas=1,nBas2
     trace3=trace3+Chi_ao_iw_v(ibas,ibas)
    enddo
   
    ! For RPA
    call diagonalize_general_matrix(nBas2,eps,eigval_eps,eigenv_eps)
    do ibas=1,nBas2
     trace2=trace2+Log(abs(eigval_eps(ibas))) ! Recall that Tr [ U Log[Eig_Mat] U^T ] = Tr [ U^T U Log[Eig_Mat] ] =  Tr [ Log[Eig_Mat] ]
                                              ! (i.e., we just need to sum the Log of the eigenvalues)
    enddo
   
    ! Compute EcRPA and EcGM from traces 
    EcRPA=EcRPA+wweight(ifreq)*(trace2+trace)/(2d0*pi)
    EcGM =EcGM -wweight(ifreq)*(trace3-trace)/(2d0*pi)
   
    ! Building Wp in AO basis
    Wp_ao_iw(:,:)=matmul(vMat(:,:),Chi_ao_iw_v(:,:))
    if(ifreq==1 .and. fulltest) then
      do ibas=1,nBas
       do jbas=1,nBas
        do kbas=1,nBas
         do lbas=1,nBas
          Wp_AO(ibas,jbas,kbas,lbas)=Wp_ao_iw(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)
         enddo
        enddo
       enddo
      enddo
      call AOtoMO_ERI_RHF(nBas,nOrb,cHFB,Wp_AO,Wp_MO)
      do iorb=1,nOrb
       do jorb=1,nOrb
        do korb=1,nOrb
         do lorb=1,nOrb
          Wp_tmp(1+(korb-1)+(iorb-1)*nBas,1+(lorb-1)+(jorb-1)*nBas)=Wp_MO(iorb,jorb,korb,lorb)
         enddo
        enddo
       enddo
      enddo
      write(*,*) ' ' 
      write(*,*) 'HFB Wp_MO(i w1) ' 
      write(*,*) ' ' 
      do iorb=1,nOrb*nOrb
       write(*,'(*(f10.5))') Wp_tmp(iorb,:)
      enddo
      write(*,*) ' ' 
    endif

    ! Build G(iw+wtest)
    eta=0d0
    ! Ghe
     Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
     Mat2(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
     Mat3(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
     Mat4(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
     weval=wtest+im*wcoord(ifreq)
     call G_AO_HFB(nBas,nOrb,nOrb_twice,eta,cHFB,eHFB,weval,Mat1,Mat2,Mat3,Mat4,G_ao_1)
     weval=wtest-im*wcoord(ifreq)
     call G_AO_HFB(nBas,nOrb,nOrb_twice,eta,cHFB,eHFB,weval,Mat1,Mat2,Mat3,Mat4,G_ao_2)
    ! Ghh
     Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
     Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
     Mat3(1:nOrb,1:nOrb)=-U_QP(nOrb+1:nOrb_twice,1:nOrb)
     Mat4(1:nOrb,1:nOrb)= U_QP(1:nOrb,1:nOrb)
     weval=wtest+im*wcoord(ifreq)
     call G_AO_HFB(nBas,nOrb,nOrb_twice,eta,cHFB,eHFB,weval,Mat1,Mat2,Mat3,Mat4,G_ao_3)
     weval=wtest-im*wcoord(ifreq)
     call G_AO_HFB(nBas,nOrb,nOrb_twice,eta,cHFB,eHFB,weval,Mat1,Mat2,Mat3,Mat4,G_ao_4)

    ! Sigma_he/hh_c(wtest)
    do ibas=1,nBas
     do jbas=1,nBas
      do kbas=1,nBas
       do lbas=1,nBas
        ! This is  G^he W
        Sigma_he_c_ao(ibas,jbas)=Sigma_he_c_ao(ibas,jbas)-(G_ao_1(kbas,lbas)+G_ao_2(kbas,lbas))  &
                                *Wp_ao_iw(1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)     &
                                *wweight(ifreq)/(2d0*pi)
        ! This is -G^hh W
        Sigma_hh_c_ao(ibas,jbas)=Sigma_hh_c_ao(ibas,jbas)+(G_ao_3(kbas,lbas)+G_ao_4(kbas,lbas))  &
                                *Wp_ao_iw(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)     &
                                *wweight(ifreq)/(2d0*pi)
       enddo
      enddo
     enddo
    enddo
   
  enddo

  ! Print Sigma_he/hh_c_ao
  if(fulltest) then

   write(*,*) ' '
   write(*,'(a,f15.8,a,f15.8,a)') ' HFB Sigma_he_c(wtest) in AO for wtest=(',Real(wtest),",",Aimag(wtest),")"
   write(*,*) ' '
   do iorb=1,nOrb
    write(*,'(*(f10.5))') Real(Sigma_he_c_ao(iorb,:))
   enddo
   write(*,*) ' '
   write(*,'(a,f15.8,a,f15.8,a)') ' HFB Sigma_hh_c(wtest) in AO for wtest=(',Real(wtest),",",Aimag(wtest),")"
   write(*,*) ' '
   do iorb=1,nOrb
    write(*,'(*(f10.5))') Real(Sigma_hh_c_ao(iorb,:))
   enddo
   write(*,*) ' '

  endif

  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A60,F15.6,A3)') '         phRPA correlation energy = ',EcRPA,' au'
  write(*,'(2X,A60,F15.6,A3)') '            GM correlation energy = ',EcGM,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  ! Deallocate arrays
  deallocate(eps,epsm1,eigval_eps,eigenv_eps)
  deallocate(Chi0_ao_iw)
  deallocate(Chi0_ao_iw_v,Chi_ao_iw_v,Wp_ao_iw)
  deallocate(Wp_AO,Wp_MO,Wp_tmp)
  deallocate(vMat)
  deallocate(Sigma_he_c_ao,G_ao_1,G_ao_2)
  deallocate(Sigma_hh_c_ao,G_ao_3,G_ao_4)
  deallocate(Mat1,Mat2,Mat3,Mat4)

  call wall_time(end_Xoiw)
  t_Xoiw = end_Xoiw - start_Xoiw
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Xo(i w) test = ',t_Xoiw,' seconds'
  write(*,*)

end subroutine

