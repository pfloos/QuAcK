subroutine Gitau2Chi0iw_mo_HFB(nOrb,nOrb_twice,verbose,eHFB,ntimes,wcoord,U_QP,Chi0_mo_iw)

! Restricted Xo(i tau) [ and Xo(i w) ] computed from G(i tau)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: ntimes
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: wcoord
  double precision,intent(in)   :: eHFB(nOrb_twice)
  double precision,intent(in)   :: U_QP(nOrb_twice,nOrb_twice)

! Local variables

  logical                       :: lesser

  integer                       :: kind_int,itau
  integer                       :: iorb,jorb,korb,lorb

  double precision              :: start_Gitau2Xoiw     ,end_Gitau2Xoiw       ,t_Gitau2Xoiw

  double precision              :: alpha,beta,lim_inf,lim_sup
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: Mat1(:,:),Mat2(:,:)

  complex*16                    :: product1
  complex*16                    :: product2
  complex*16,allocatable        :: Glesser_he(:,:),Ggreater_he(:,:)
  complex*16,allocatable        :: Glesser_hh(:,:),Ggreater_ee(:,:)
  complex*16,allocatable        :: Chi0_mo_itau(:,:)

! Output variables
  complex*16,intent(out)        :: Chi0_mo_iw(nOrb*nOrb,nOrb*nOrb)
  
!------------------------------------------------------------------------
! Build G(i tau) in MO orbis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------

 call wall_time(start_Gitau2Xoiw)

 if(verbose/=0) then
  write(*,*)     
  write(*,*)'*********************************************'
  write(*,*)'* HFB G(i tau) and Xo(i w) construction     *'
  write(*,*)'*********************************************'
  write(*,*)
 endif
 
 Chi0_mo_iw(:,:)=czero

 allocate(Mat1(nOrb,nOrb),Mat2(nOrb,nOrb)) 
 allocate(Glesser_he(nOrb,nOrb),Ggreater_he(nOrb,nOrb)) 
 allocate(Glesser_hh(nOrb,nOrb),Ggreater_ee(nOrb,nOrb)) 
 allocate(Chi0_mo_itau(nOrb*nOrb,nOrb*nOrb)) 

 Chi0_mo_itau=czero

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

  ! Ghe
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  lesser=.true.  ! Build MO G<
  call build_Glorg_mo_HFB(nOrb,nOrb_twice,tcoord(itau),Glesser_he,eHFB,Mat1,Mat2,lesser)
  Mat1(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
  lesser=.false. ! Build MO G>
  call build_Glorg_mo_HFB(nOrb,nOrb_twice,tcoord(itau),Ggreater_he,eHFB,Mat1,Mat2,lesser)
  ! Ghh and Gee
  Mat1(1:nOrb,1:nOrb)=U_QP(1:nOrb,1:nOrb)
  Mat2(1:nOrb,1:nOrb)=U_QP(nOrb+1:nOrb_twice,1:nOrb)
  lesser=.true.  ! Build MO G<
  call build_Glorg_mo_HFB(nOrb,nOrb_twice,tcoord(itau),Glesser_hh,eHFB,Mat1,Mat2,lesser)
  lesser=.false. ! Build MO G> [this one should be computed with -Mat2 according to Eq. 45 in V. Soma et al Phys. Rev. C, 84, 064317 (2011)] 
  call build_Glorg_mo_HFB(nOrb,nOrb_twice,tcoord(itau),Ggreater_ee,eHFB,Mat1,Mat2,lesser)

  ! Xo(i tau) = -2i G<(i tau) G>(-i tau) [ We store the spinless elements in this manner ]
  do iorb=1,nOrb
   do jorb=1,nOrb
    do korb=jorb,nOrb
     do lorb=iorb,nOrb                       
                          ! r1   r2'               r2   r1'
      product1 = Glesser_he(iorb,jorb)*Ggreater_he(korb,lorb) &
               - Glesser_hh(iorb,jorb)*Ggreater_ee(korb,lorb)   ! [we add a minus to compensate not using -Mat2 according to Soma's paper]
      product2 = Glesser_he(korb,lorb)*Ggreater_he(iorb,jorb) &
               - Glesser_hh(korb,lorb)*Ggreater_ee(iorb,jorb)   ! [we add a minus to compensate not using -Mat2 according to Soma's paper]
      if(abs(product1)<1e-12) product1=czero
      if(abs(product2)<1e-12) product2=czero
      Chi0_mo_itau(1+(lorb-1)+(iorb-1)*nOrb,1+(korb-1)+(jorb-1)*nOrb) = product1+product2
     enddo
    enddo
   enddo
  enddo
  Chi0_mo_itau=-2d0*im*Chi0_mo_itau ! The 2 factor is added to account for both spin contributions 
                                    ! [ i.e., for Ghe Ghe take (up,up,up,up) and (down,down,down,down) 
                                    !                                 & 
                                    !       for Ghh Gee (up,down,down,up) and (down,up,up,down) ]

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

 ! Deallocate dyn arrays
 deallocate(tcoord,tweight) 
 deallocate(Mat1,Mat2) 
 deallocate(Glesser_he,Ggreater_he) 
 deallocate(Glesser_hh,Ggreater_ee) 
 deallocate(Chi0_mo_itau) 

end subroutine 


subroutine build_Glorg_mo_HFB(nOrb,nOrb_twice,tau,Glorg_mo,eHFB,Mat1,Mat2,lesser)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: lesser

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nOrb_twice

  double precision,intent(in)   :: tau
  double precision,intent(in)   :: eHFB(nOrb_twice)
  double precision,intent(in)   :: Mat1(nOrb,nOrb)
  double precision,intent(in)   :: Mat2(nOrb,nOrb)

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
  
  do iorb=1,nOrb
    Glorg_mo(:,:) = Glorg_mo(:,:) + fact*im*Exp(eHFB(iorb)*tau) &
              * matmul(Mat1(:,iorb:iorb),transpose(Mat2(:,iorb:iorb)))
  enddo

end subroutine

