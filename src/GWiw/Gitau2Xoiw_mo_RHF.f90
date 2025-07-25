subroutine Gitau2Chi0iw_mo_RHF(nOrb,nO,verbose,eHF,ntimes,wcoord,Chi0_mo_iw)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: wcoord

! Local variables

  logical                       :: lesser

  integer                       :: kind_int,itau
  integer                       :: iorb,jorb,korb,lorb

  double precision              :: start_Gitau2Xoiw     ,end_Gitau2Xoiw       ,t_Gitau2Xoiw

  double precision              :: chem_pot
  double precision              :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable  :: tweight(:),tcoord(:)

  complex*16                    :: product
  complex*16,allocatable        :: Glesser(:,:),Ggreater(:,:)
  complex*16,allocatable        :: Chi0_mo_itau(:,:)

! Output variables
  double precision,intent(inout):: eHF(nOrb)
  complex*16,intent(out)        :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)
  
!------------------------------------------------------------------------
! Build G(i tau) in MO orbis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------
 
 call wall_time(start_Gitau2Xoiw)

 if(verbose/=0) then
  write(*,*)     
  write(*,*)'*********************************************'
  write(*,*)'* RHF G(i tau) and Xo(i w) construction     *'
  write(*,*)'*********************************************'
  write(*,*)
 endif

 chem_pot = 0.5d0*(eHF(nO)+eHF(nO+1))
 if(verbose/=0) then
  write(*,'(A33,1X,F16.10,A3)') ' Chemical potential  = ',chem_pot,' au'
 endif
 eHF(:) = eHF(:)-chem_pot
 Chi0_mo_iw(:,:)=czero
   
 allocate(Glesser(nOrb,nOrb),Ggreater(nOrb,nOrb)) 
 allocate(Chi0_mo_itau(nOrb*nOrb,nOrb*nOrb)) 
 Chi0_mo_itau(:,:)=czero

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
 if(verbose/=0) then
  write(*,*)
  write(*,*) '    ----------------------------'
  write(*,'(A28,1X)') 'Testing the time quadrature'
  write(*,'(A28, I12)') 'Number of times (grid)= ', ntimes
  write(*,'(A28,1X,F16.10)') 'PI value error',abs(alpha-acos(-1d0))
  write(*,*) '    ----------------------------'
  write(*,*)
 endif
  alpha = 0d0; beta = 0d0;

 ! time grid Xo(i tau) = -2i G<(i tau) G>(-i tau)
 !  and Fourier transform Xo(i tau) -> Xo(i w)
 do itau=1,ntimes

  lesser=.true.  ! Build MO G<
  call build_Glorg_mo_RHF(nOrb,nO,tcoord(itau),Glesser,eHF,lesser)
  lesser=.false. ! Build MO G>
  call build_Glorg_mo_RHF(nOrb,nO,tcoord(itau),Ggreater,eHF,lesser)

  ! Xo(i tau) = -2i G<(i tau) G>(-i tau) [ In HF using the MO basis, Ggreater and Glesser are diagonal ]
  do iorb=1,nO
   jorb=iorb
   do korb=nO+1,nOrb
    lorb=korb                       
                     ! r1   r2'            r2   r1'
    product = Glesser(iorb,jorb)*Ggreater(korb,lorb)
    if(abs(product)<1e-12) product=czero
    Chi0_mo_itau(1+(lorb-1)+(iorb-1)*nOrb,1+(korb-1)+(jorb-1)*nOrb) = product
   enddo
  enddo
  Chi0_mo_itau=-2d0*im*Chi0_mo_itau ! The 2 factor is added to account for both spin contributions [ i.e., (up,up,up,up) and (down,down,down,down) ]

  ! Xo(i tau) -> Xo(i w)
  Chi0_mo_iw(:,:) = Chi0_mo_iw(:,:) - im*tweight(itau)*Chi0_mo_itau(:,:)*Exp(im*tcoord(itau)*wcoord)

 enddo

 ! Complete the Xo(i tau) -> Xo(i w)
 Chi0_mo_iw(:,:) = 2d0*Real(Chi0_mo_iw(:,:))

 call wall_time(end_Gitau2Xoiw)
 
 t_Gitau2Xoiw = end_Gitau2Xoiw - start_Gitau2Xoiw
 if(verbose/=0) then
  write(*,'(A65,1X,F9.3,A8)') 'Total wall time for Gitau2Chi0iw = ',t_Gitau2Xoiw,' seconds'
  write(*,*)
 endif

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot
 deallocate(Glesser,Ggreater,tcoord,tweight) 
 deallocate(Chi0_mo_itau) 

end subroutine 


subroutine build_Glorg_mo_RHF(nOrb,nO,tau,Glorg_mo,eHF,lesser)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: lesser

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  integer                       :: iorb

  double precision              :: chem_pot,fact

! Output variables
  complex*16,intent(out)        :: Glorg_mo(nOrb,nOrb)
  
!------------------------------------------------------------------------
! Build G<(i tau) and G>(i tau)
!------------------------------------------------------------------------

  if(lesser) then ! G<
   fact=1d0
  else            ! G>
   fact=-1d0
  endif
  Glorg_mo=czero
  
  if(lesser) then ! G<
   do iorb=1,nO
     Glorg_mo(iorb,iorb) = fact*im*Exp(eHF(iorb)*tau)
   enddo
  else            ! G>
   do iorb=nO+1,nOrb
     Glorg_mo(iorb,iorb) = fact*im*Exp(-eHF(iorb)*tau)
   enddo
  endif 

end subroutine

