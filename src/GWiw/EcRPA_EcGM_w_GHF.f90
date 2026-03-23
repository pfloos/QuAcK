subroutine EcRPA_EcGM_w_GHF(nBas,nBas2,nO,verbose,e_GHF,cGHF,nfreqs,wweight,wcoord,vMAT,EGHF)

! Restricted EcRPA and EcGM for GHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nO

  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: cGHF(nBas2,nBas2)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)

! Local variables

  integer                       :: ntimes
  integer                       :: ntimes_twice
  integer                       :: ifreq,itau
  integer                       :: ibas,jbas,kbas,lbas
  integer                       :: qbas,rbas,sbas,tbas
  integer                       :: nBasSq

  double precision              :: eta
  double precision              :: chem_pot
  double precision              :: EcRPA,EcGM
  double precision              :: trace1,trace2,trace3
  double precision,allocatable  :: Tmp_ao_w(:,:)
  double precision,allocatable  :: eigval_Xov(:)
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: sint2w_weight(:,:)
  double precision,allocatable  :: cost2w_weight(:,:)
  double precision,allocatable  :: cosw2t_weight(:,:)
  double precision,allocatable  :: sinw2t_weight(:,:)

  complex*16                    :: product
  complex*16,allocatable        :: Chi0_ao_itau(:,:)
  complex*16,allocatable        :: Chi0_ao_iw(:,:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_minus_itau(:,:),G_plus_itau(:,:)

! Ouput variables

  double precision,intent(inout):: e_GHF(nBas2)

!

  ntimes=nfreqs
  ntimes_twice=2*ntimes
  nBasSq=nBas*nBas
  EcRPA=0d0; EcGM=0d0;
  allocate(eigval_Xov(nBasSq))
  allocate(G_ao_itau(ntimes_twice,nBas2,nBas2))
  allocate(G_minus_itau(nBas2,nBas2),G_plus_itau(nBas2,nBas2)) 
  allocate(Chi0_ao_iw(nfreqs,nBasSq,nBasSq),Tmp_ao_w(nBasSq,nBasSq))
  allocate(Chi0_ao_itau(nBasSq,nBasSq))
  allocate(tweight(ntimes),tcoord(ntimes))
  allocate(sint2w_weight(nfreqs,ntimes))
  allocate(cost2w_weight(nfreqs,ntimes))
  allocate(cosw2t_weight(ntimes,nfreqs))
  allocate(sinw2t_weight(ntimes,nfreqs))
  chem_pot = 0.5d0*(e_GHF(nO)+e_GHF(nO+1))
  e_GHF(:)=e_GHF(:)-chem_pot

 ! Read grids 
  call read_scGX_grids(ntimes,nfreqs,tcoord,tweight,wcoord,wweight,sint2w_weight,cost2w_weight, &
                       cosw2t_weight,sinw2t_weight,0)

 ! Build Go(i tau)
  do itau=1,ntimes
   call G_AO_GHF_itau(nBas2,nO, tcoord(itau),G_plus_itau ,cGHF,e_GHF)
   call G_AO_GHF_itau(nBas2,nO,-tcoord(itau),G_minus_itau,cGHF,e_GHF)
   G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:)
   G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)
  enddo

 ! Build using the time grid Xo(i tau) = -i \sum_ss' Gss'(i tau) Gs's(-i tau)
 !  then Fourier transform Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:)=czero
  do itau=1,ntimes
   ! Xo(i tau) = -i \sum_ss' Gss'(i tau) Gs's(-i tau)
   do ibas=1,nBas
    qbas=nBas+1+(ibas-1)
    do jbas=1,nBas
     rbas=nBas+1+(jbas-1)
     do kbas=1,nBas
      sbas=nBas+1+(kbas-1)
      do lbas=1,nBas
       tbas=nBas+1+(lbas-1)
                                   ! r1   r2'                    r2   r1'
       product = G_ao_itau(2*itau-1,ibas,jbas)*G_ao_itau(2*itau,kbas,lbas) &
               + G_ao_itau(2*itau-1,ibas,rbas)*G_ao_itau(2*itau,sbas,lbas) &
               + G_ao_itau(2*itau-1,qbas,jbas)*G_ao_itau(2*itau,kbas,tbas) &
               + G_ao_itau(2*itau-1,qbas,rbas)*G_ao_itau(2*itau,sbas,tbas)
       if(abs(product)<1e-12) product=czero
       Chi0_ao_itau(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
      enddo
     enddo
    enddo
   enddo
   Chi0_ao_itau=-im*Chi0_ao_itau 
   ! Xo(i tau) -> Xo(i w) [ the weight already contains the cos(tau w) and a factor 2 because int_-Infty ^Infty -> 2 int_0 ^Infty ]
   do ifreq=1,nfreqs
    Chi0_ao_iw(ifreq,:,:) = Chi0_ao_iw(ifreq,:,:) - im*cost2w_weight(ifreq,itau)*Chi0_ao_itau(:,:)
   enddo
  enddo
  ! Complete the Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:) = Real(Chi0_ao_iw(:,:,:)) ! The factor 2 is stored in the weight [ and we just retain the real part ]


! Imaginary freqs contribution

  do ifreq=1,nfreqs
   
   ! Initialize
   trace1=0d0; trace2=0d0; trace3=0d0;
   Tmp_ao_w=0d0

   ! Tr [ Xo v ] and define eps
   Tmp_ao_w(:,:)=matmul(Real(Chi0_ao_iw(ifreq,:,:)),vMAT(:,:))
   do ibas=1,nBasSq
    trace1=trace1+Tmp_ao_w(ibas,ibas)
   enddo

   ! Diagonalize Xo v
   call diagonalize_general_matrix(nBasSq,Tmp_ao_w,eigval_Xov,Tmp_ao_w)

   ! GM 
   ! Tr [ X v ] = Tr [ (1 - Xo v)^-1 Xo v ] = Tr [ (I - U^-1 Xo v U)^-1 (U^-1 Xo v U) ] = Tr [ ( I - Eigval_mat)^-1 Eigval_mat ]
   !            = sum_p eigval_p / (1-eigval_p)
   ! RPA 
   ! Tr [ Log(I - Xo v) ] = Tr [ U^-1 Log(I - U^-1 Xo v U) U ] = Tr [ U U^-1 Log(I- Eigval_mat) ] = Tr [ Log(I - Eigval_mat) ]
   !            = sum_p Log(1-eigval_p)
   do ibas=1,nBasSq
    trace2=trace2+Log(abs(1d0-eigval_Xov(ibas)))
    trace3=trace3+eigval_Xov(ibas)/(1d0-eigval_Xov(ibas))
   enddo

   ! Compute EcRPA and EcGM from traces 
   EcRPA=EcRPA+wweight(ifreq)*(trace2+trace1)/(2d0*pi)
   EcGM =EcGM -wweight(ifreq)*(trace3-trace1)/(2d0*pi)

  enddo
  e_GHF(:)=e_GHF(:)+chem_pot

! Print results
 
  if(verbose/=0) then 
   write(*,*)
   write(*,*) '******************************************************************'
   write(*,*) '* Generalized EcRPA and EcGM computed with imaginary frequencies *'
   write(*,*) '******************************************************************'
   write(*,*)
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,'(2X,A60,F15.6,A3)') '       G-phRPA correlation energy = ',EcRPA,' au'
   write(*,'(2X,A60,F15.6,A3)') '       G-phRPA total energy       = ',EGHF+EcRPA,' au'
   write(*,'(2X,A60,F15.6,A3)') '          G-GM correlation energy = ',EcGM,' au'
   write(*,'(2X,A60,F15.6,A3)') '          G-GM total energy       = ',EGHF+EcGM,' au'
   write(*,*)'-------------------------------------------------------------------------------'
   write(*,*)
  endif

  ! Deallocate arrays
  deallocate(tcoord,tweight) 
  deallocate(sint2w_weight)
  deallocate(cost2w_weight)
  deallocate(cosw2t_weight)
  deallocate(sinw2t_weight)
  deallocate(G_ao_itau)
  deallocate(G_minus_itau,G_plus_itau) 
  deallocate(Chi0_ao_iw,Chi0_ao_itau)
  deallocate(Tmp_ao_w)
  deallocate(eigval_Xov)

end subroutine

