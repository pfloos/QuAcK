subroutine iGtau2Chi0iw_HFB(nBas,nBas2,nOrb,nOrb2,cHFB,eHFB,nfreqs,ntimes,wcoord,U_QP,Chi0_ao_iw)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb2

  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: eHFB(nOrb2)
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: U_QP(nOrb2,nOrb2)

! Local variables

  logical                       :: lesser

  integer                       :: kind_int,itau,ifreq
  integer                       :: ibas,jbas,kbas,lbas

  double precision              :: start_iG     ,end_iG       ,t_iG

  double precision              :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: Mat1(:,:),Mat2(:,:)

  complex*16                    :: product
  complex*16,allocatable        :: Glesser_he(:,:),Ggreater_he(:,:)
  complex*16,allocatable        :: Glesser_hh(:,:),Ggreater_ee(:,:)
  complex*16,allocatable        :: Chi0_ao_itau(:,:)

! Output variables
  complex*16,intent(out)        :: Chi0_ao_iw(nfreqs,nBas*nBas,nBas*nBas)
  
!------------------------------------------------------------------------
! Build iG(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------

 write(*,*)     
 write(*,*)'*****************************************'
 write(*,*)'* G(i tau) and Xo(i w) construction     *'
 write(*,*)'*****************************************'
 write(*,*)
 
 call wall_time(start_iG)

 Chi0_ao_iw(:,:,:)=czero
   
 allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb)) 
 allocate(Glesser_he(nBas,nBas),Ggreater_he(nBas,nBas)) 
 allocate(Glesser_hh(nBas,nBas),Ggreater_ee(nBas,nBas)) 
 allocate(Chi0_ao_itau(nBas2,nBas2)) 

!-------------------------!
! Prepare time Quadrature !
!-------------------------!
  kind_int = 1
  lim_inf = 0d0; lim_sup = 1d0;
  alpha = 0d0;  beta  = 0d0;
  allocate(tweight(ntimes),tcoord(ntimes))
  call cgqf(ntimes,kind_int,alpha,beta,lim_inf,lim_sup,tcoord,tweight)
  tweight(:)=tweight(:)/((1d0-tcoord(:))**2d0)
  tcoord(:)=tcoord(:)/(1d0-tcoord(:))
  ! Check how good we integrate for beta = eA-eI
  alpha = 0d0; beta = abs(eHFB(nOrb)+eHFB(1));
  do itau=1,ntimes
   alpha=alpha+tweight(itau)*(beta*2d0/(beta**2d0+tcoord(itau)**2d0))
  enddo
  write(*,*)
  write(*,*) '    ----------------------------'
  write(*,'(A28,1X)') 'Testing the time quadrature'
  write(*,'(A28, I12)') 'Number of times (grid)= ', ntimes
  write(*,'(A28,1X,F16.10)') 'PI value error',abs(alpha-acos(-1d0))
  write(*,*) '    ----------------------------'
  write(*,*)
  alpha = 0d0; beta = 0d0;

 ! time grid Xo(i tau) = -2i G<(i tau) G>(-i tau)
 !  and Fourier transform Xo(i tau) -> Xo(i w)
 do itau=1,ntimes

  ! Ghe
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  lesser=.true.  ! Build AO G<
  call build_Glorg_HFB(nBas,nOrb,nOrb2,tcoord(itau),Glesser_he,cHFB,eHFB,Mat1,Mat2,lesser)
  Mat1(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb2,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb2,1:nOrb)
  lesser=.false. ! Build AO G>
  call build_Glorg_HFB(nBas,nOrb,nOrb2,tcoord(itau),Ggreater_he,cHFB,eHFB,Mat1,Mat2,lesser)
  ! Ghh and Gee
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb2,1:nOrb)
  lesser=.true.  ! Build AO G<
  call build_Glorg_HFB(nBas,nOrb,nOrb2,tcoord(itau),Glesser_hh,cHFB,eHFB,Mat1,Mat2,lesser)
  lesser=.false. ! Build AO G>
  call build_Glorg_HFB(nBas,nOrb,nOrb2,tcoord(itau),Ggreater_ee,cHFB,eHFB,Mat1,Mat2,lesser)

  ! Xo(i tau) = -2i G<(i tau) G>(-i tau)
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas                       
                          ! r1   r2'               r2   r1'
      product = Glesser_he(ibas,jbas)*Ggreater_he(kbas,lbas) &
              - Glesser_hh(ibas,jbas)*Ggreater_ee(kbas,lbas)
      if(abs(product)<1e-12) product=czero
      Chi0_ao_itau(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
     enddo
    enddo
   enddo
  enddo
  Chi0_ao_itau=-2d0*im*Chi0_ao_itau ! The 2 factor is added to account for both spin contributions [ i.e., (up,up) and (down,down) ]

  ! Xo(i tau) -> Xo(i w)
  do ifreq=1,nfreqs
    Chi0_ao_iw(ifreq,:,:) = Chi0_ao_iw(ifreq,:,:) - im*tweight(itau)*Chi0_ao_itau(:,:)*Exp(im*tcoord(itau)*wcoord(ifreq))
  enddo 

 enddo

 ! Complete the Xo(i tau) -> Xo(i w)
 Chi0_ao_iw(:,:,:) = 2d0*Real(Chi0_ao_iw(:,:,:))

 call wall_time(end_iG)
 
 t_iG = end_iG - start_iG
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for iGtau = ',t_iG,' seconds'
 write(*,*)

 ! Deallocate dyn arrays
 deallocate(tcoord,tweight) 
 deallocate(Mat1,Mat2) 
 deallocate(Glesser_he,Ggreater_he) 
 deallocate(Glesser_hh,Ggreater_ee) 
 deallocate(Chi0_ao_itau) 

end subroutine 


subroutine build_Glorg_HFB(nBas,nOrb,nOrb2,tau,Glorg,cHFB,eHFB,Mat1,Mat2,lesser)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: lesser

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb2

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: cHFB(nBas,nOrb)
  double precision,intent(in)   :: eHFB(nOrb2)
  double precision,intent(in)   :: Mat1(nOrb,nOrb)
  double precision,intent(in)   :: Mat2(nOrb,nOrb)

! Local variables

  double precision              :: start_iG     ,end_iG       ,t_iG

  integer                       :: iorb

  double precision              :: chem_pot,fact

  complex*16,allocatable        :: Gtmp(:,:)

! Output variables
  complex*16,intent(out)        :: Glorg(nBas,nBas)
  
!------------------------------------------------------------------------
! Build G<(i tau) and G>(i tau)
!------------------------------------------------------------------------

  if(lesser) then ! G<
   fact=1d0
  else            ! G>
   fact=-1d0
  endif
  allocate(Gtmp(nOrb,nOrb)) 
  Gtmp=czero
  
  if(lesser) then ! G<
   do iorb=1,nOrb
     Gtmp(:,:) = Gtmp(:,:) + fact*im*Exp(eHFB(iorb)*tau) &
               * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
   enddo
  else            ! G>
   do iorb=1,nOrb
     Gtmp(:,:) = Gtmp(:,:) + fact*im*Exp(eHFB(iorb)*tau) &
               * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
   enddo
  endif 

  Glorg=matmul(matmul(cHFB,Gtmp),transpose(cHFB))
  
  deallocate(Gtmp)

end subroutine

