subroutine Gitau2Chi0iw_RHF(nBas,nOrb,nO,cHF,eHF,nfreqs,ntimes,wcoord,Chi0_ao_iw)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: cHF(nBas,nOrb)

! Local variables

  logical                       :: lesser

  integer                       :: kind_int,itau,ifreq
  integer                       :: ibas,jbas,kbas,lbas,nBas2

  double precision              :: start_Gitau2Xoiw     ,end_Gitau2Xoiw       ,t_Gitau2Xoiw

  double precision              :: chem_pot
  double precision              :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable  :: tweight(:),tcoord(:)

  complex*16                    :: product
  complex*16,allocatable        :: Glesser(:,:),Ggreater(:,:)
  complex*16,allocatable        :: Chi0_ao_itau(:,:)

! Output variables
  double precision,intent(inout):: eHF(nOrb)
  complex*16,intent(out)        :: Chi0_ao_iw(nfreqs,nBas*nBas,nBas*nBas)
  
!------------------------------------------------------------------------
! Build G(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------

 write(*,*)     
 write(*,*)'*********************************************'
 write(*,*)'* RHF G(i tau) and Xo(i w) construction     *'
 write(*,*)'*********************************************'
 write(*,*)
 
 call wall_time(start_Gitau2Xoiw)

 nBas2=nBas*nBas
 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
 eHF(:) = eHF(:)-chem_pot
 Chi0_ao_iw(:,:,:)=czero
   
 allocate(Glesser(nBas,nBas),Ggreater(nBas,nBas)) 
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
  alpha = 0d0; beta = eHF(nOrb)-eHF(1);
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

  lesser=.true.  ! Build AO G<
  call build_Glorg_RHF(nBas,nOrb,nO,tcoord(itau),Glesser,cHF,eHF,lesser)
  lesser=.false. ! Build AO G>
  call build_Glorg_RHF(nBas,nOrb,nO,tcoord(itau),Ggreater,cHF,eHF,lesser)

  ! Xo(i tau) = -2i G<(i tau) G>(-i tau)
  do ibas=1,nBas
   do jbas=1,nBas
    do kbas=1,nBas
     do lbas=1,nBas                       
                       ! r1   r2'            r2   r1'
      product = Glesser(ibas,jbas)*Ggreater(kbas,lbas)
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

 call wall_time(end_Gitau2Xoiw)
 
 t_Gitau2Xoiw = end_Gitau2Xoiw - start_Gitau2Xoiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Gitau2Chi0iw = ',t_Gitau2Xoiw,' seconds'
 write(*,*)

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot
 deallocate(Glesser,Ggreater,tcoord,tweight) 
 deallocate(Chi0_ao_itau) 

end subroutine 


subroutine build_Glorg_RHF(nBas,nOrb,nO,tau,Glorg,cHF,eHF,lesser)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: lesser

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

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
   do iorb=1,nO
     Gtmp(iorb,iorb) = fact*im*Exp(eHF(iorb)*tau)
   enddo
  else            ! G>
   do iorb=nO+1,nOrb
     Gtmp(iorb,iorb) = fact*im*Exp(-eHF(iorb)*tau)
   enddo
  endif 

  Glorg=matmul(matmul(cHF,Gtmp),transpose(cHF))
  
  deallocate(Gtmp)

end subroutine

