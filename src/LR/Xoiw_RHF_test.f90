subroutine build_Xoiw_RHF_test(nBas,nOrb,nO,cHF,eHF,nfreqs,ntimes,wweight,wcoord,ERI_AO)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ibas,jbas,kbas,lbas,ifreq
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: nBas2

  double precision              :: start_Xoiw   ,end_Xoiw     ,t_Xoiw
  double precision              :: EcRPA,EcGM,trace,trace2,trace3
  double precision,allocatable  :: Chi0_ao_iw_v(:,:)
  double precision,allocatable  :: Chi_ao_iw_v(:,:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: eigval_eps(:)
  double precision,allocatable  :: eps(:,:)
  double precision,allocatable  :: epsm1(:,:)
  double precision,allocatable  :: eigenv_eps(:,:)
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: Wp_tmp(:,:)
  double precision,allocatable  :: Wp_AO(:,:,:,:)
  double precision,allocatable  :: Wp_MO(:,:,:,:)

  complex *16,allocatable       :: Chi0_ao_iw(:,:,:)
!

  nBas2=nBas*nBas

!------------------------------------------------------------------------
! Build G(i tau) in AO basis
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'*******************************************'
  write(*,*)'* Use RHF Xo(i w) to build RHF X(i w),    *'
  write(*,*)'*       compute EcRPA and EcGM,           *'
  write(*,*)'*          and build Wp (i w)             *'
  write(*,*)'*******************************************'
  write(*,*)

  allocate(eps(nBas2,nBas2),epsm1(nBas2,nBas2),eigenv_eps(nBas2,nBas2),eigval_eps(nBas2))
  allocate(Chi0_ao_iw(nfreqs,nBas2,nBas2))
  allocate(Chi0_ao_iw_v(nBas2,nBas2),Chi_ao_iw_v(nBas2,nBas2))
  allocate(Wp_AO(nBas,nBas,nBas,nBas),Wp_MO(nOrb,nOrb,nOrb,nOrb),Wp_tmp(nOrb*nOrb,nOrb*nOrb))
  allocate(Wp_ao_iw(nBas2,nBas2))
  allocate(vMAT(nBas2,nBas2))

!-----------------------------!
! Store v also as a 2D matrix !
!-----------------------------!
 
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas
      vMAT(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)=ERI_AO(ibas,jbas,kbas,lbas)
     enddo
    enddo
   enddo
  enddo

!-----------------------------!
! Build Xo(i w) from G(i tau) !
!-----------------------------!

  call Gitau2Chi0iw_RHF(nBas,nOrb,nO,cHF,eHF,nfreqs,ntimes,wcoord,Chi0_ao_iw)

  call wall_time(start_Xoiw)

!----------------------!
! Use Xo(i w) as usual !
!----------------------!

  EcRPA=0d0; EcGM=0d0;
  do ifreq=1,nfreqs

    ! Initialization
    trace=0d0; trace2=0d0; trace3=0d0;
    eigval_eps(:)=0d0
    eigenv_eps(:,:)=0d0
    eps(:,:)=0d0
   
    ! Build Xo v
    Chi0_ao_iw_v(:,:)=matmul(Real(Chi0_ao_iw(ifreq,:,:)),vMAT(:,:))
   
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
    Chi_ao_iw_v(:,:)=matmul(Chi_ao_iw_v(:,:),vMAT(:,:))
   
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
    Wp_ao_iw(:,:)=matmul(vMAT(:,:),Chi_ao_iw_v(:,:))
    if(ifreq==1 .and. .false.) then ! Set to true for testing ! TODO
      do ibas=1,nBas
       do jbas=1,nBas
        do kbas=1,nBas
         do lbas=1,nBas
          Wp_AO(ibas,jbas,kbas,lbas)=Wp_ao_iw(1+(kbas-1)+(ibas-1)*nBas,1+(lbas-1)+(jbas-1)*nBas)
         enddo
        enddo
       enddo
      enddo
      call AOtoMO_ERI_RHF(nBas,nOrb,cHF,Wp_AO,Wp_MO)
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
      write(*,*) 'RHF Wp_MO(i w1) ' 
      write(*,*) ' ' 
      do iorb=1,nOrb*nOrb
       write(*,'(*(f10.5))') Wp_tmp(iorb,:)
      enddo
      write(*,*) ' ' 
    endif
   
  enddo

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
  deallocate(vMAT)

  call wall_time(end_Xoiw)
  t_Xoiw = end_Xoiw - start_Xoiw
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Xo(i w) test = ',t_Xoiw,' seconds'
  write(*,*)

end subroutine

