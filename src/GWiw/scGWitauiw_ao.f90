subroutine scGWitauiw_ao(nBas,nOrb,nO,maxSCF,read_grids,ENuc,Hc,S,P_in,cHF,eHF,nfreqs,wcoord,wweight,vMAT)

! Restricted scGW

  implicit none
  include 'parameters.h'

! Input variables
 
  logical,intent(in)            :: read_grids

  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: maxSCF

  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: wweight(nfreqs)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: P_in(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: vMAT(nBas*nBas,nBas*nBas)

! Local variables
 
  logical                       :: file_exists

  integer                       :: iunit=312
  integer                       :: ntimes,ntimes_
  integer                       :: nfreqs2,nfreqs2_,ntimes_twice
  integer                       :: itau,ifreq,jfreq
  integer                       :: ibas,jbas,kbas,lbas,nBas2
  integer                       :: iter,iter_fock

  double precision              :: start_scGWitauiw     ,end_scGWitauiw       ,t_scGWitauiw

  double precision              :: Ehfl,EcGM
  double precision              :: trace1,trace2
  double precision              :: eta,diff_Pao
  double precision              :: nElectrons
  double precision              :: trace_1_rdm
  double precision              :: thrs_N,thrs_Pao
  double precision              :: chem_pot,chem_pot_saved
  double precision              :: weval(2)
  double precision,allocatable  :: tweight(:),tcoord(:)
  double precision,allocatable  :: wweight2(:),wcoord2(:)
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: Wp_ao_iw(:,:)
  double precision,allocatable  :: cHFinv(:,:)
  double precision,allocatable  :: F_ao(:,:)
  double precision,allocatable  :: P_ao(:,:)
  double precision,allocatable  :: P_ao_old(:,:)
  double precision,allocatable  :: P_ao_iter(:,:)
  double precision,allocatable  :: P_mo(:,:)
  double precision,allocatable  :: a_coef(:,:)
  double precision,allocatable  :: b_coef(:,:)

  complex*16                    :: product
  complex*16                    :: weval_cpx
  complex*16,allocatable        :: f_pq_tau(:,:)
  complex*16,allocatable        :: Sigma_c_ao(:,:,:)
  complex*16,allocatable        :: G_ao_iw2(:,:,:)
  complex*16,allocatable        :: G_ao_itau(:,:,:)
  complex*16,allocatable        :: G_ao_1(:,:),G_ao_2(:,:)
  complex*16,allocatable        :: G_minus_itau(:,:),G_plus_itau(:,:)
  complex*16,allocatable        :: Chi0_ao_itau(:,:)
  complex*16,allocatable        :: Chi0_ao_iw(:,:,:)

! Output variables
  double precision,intent(inout):: eHF(nOrb)
  
!------------------------------------------------------------------------
! Build G(i tau) in AO basis and use it to build Xo (i tau) -> Xo (i w)
!------------------------------------------------------------------------
 
 call wall_time(start_scGWitauiw)

 write(*,*)     
 write(*,*)'*******************************************'
 write(*,*)'*     scGW ( using itau and iw grids )    *'
 write(*,*)'*******************************************'
 write(*,*)

 eta=0d0
 thrs_N=1d-8
 thrs_Pao=1d-6
 nElectrons=2d0*nO
 nBas2=nBas*nBas
 chem_pot_saved = 0.5d0*(eHF(nO)+eHF(nO+1))
 chem_pot = chem_pot_saved
 Ehfl=0d0
 write(*,*)
 write(*,'(A33,1X,F16.10,A3)') ' Initial chemical potential  = ',chem_pot,' au'
 write(*,*)
 eHF(:) = eHF(:)-chem_pot_saved
   
 allocate(Chi0_ao_iw(nfreqs,nBas2,nBas2))
 allocate(a_coef(nBas,nBas),b_coef(nBas,nBas),f_pq_tau(nBas,nBas))
 allocate(P_ao(nBas,nBas),P_ao_old(nBas,nBas),P_ao_iter(nBas,nBas))
 allocate(F_ao(nBas,nBas),P_mo(nOrb,nOrb),cHFinv(nOrb,nBas),Occ(nOrb))
 allocate(G_minus_itau(nBas,nBas),G_plus_itau(nBas,nBas)) 
 allocate(G_ao_1(nBas,nBas),G_ao_2(nBas,nBas)) 
 allocate(Chi0_ao_itau(nBas2,nBas2),Wp_ao_iw(nBas2,nBas2)) 
 cHFinv=matmul(transpose(cHF),S)
 P_ao=P_in
 P_ao_iter=P_ao

!---------------!
! Prepare grids !
!---------------!

 ntimes=0;nfreqs2=nfreqs*10; ! We can use this to enforce ntimes=number1 for building a grid from [0;Infty)
 if(.not.read_grids) then
  call build_iw_itau_grid(nBas,nOrb,nO,ntimes,nfreqs2,0,cHF,eHF)
 else
  inquire(file='tcoord.txt', exist=file_exists)
  if(file_exists) then
   open(unit=iunit, form='formatted', file='tcoord.txt', status='old')
   read(iunit,*) ntimes
   close(iunit)
  endif 
  inquire(file='wcoord.txt', exist=file_exists)
  if(file_exists) then
   open(unit=iunit, form='formatted', file='wcoord.txt', status='old')
   read(iunit,*) nfreqs2
   close(iunit)
  endif 
 endif
 ntimes_twice=2*ntimes
 allocate(tweight(ntimes),tcoord(ntimes))
 allocate(wweight2(nfreqs2),wcoord2(nfreqs2))
 allocate(Sigma_c_ao(nfreqs2,nBas,nBas),G_ao_iw2(nfreqs2,nBas,nBas))
 allocate(G_ao_itau(ntimes_twice,nBas,nBas))
 write(*,*)
 write(*,'(*(a,i25))') ' Final ntimes grid size ',ntimes
 write(*,'(*(a,i25))') ' Final nfreqs grid size ',nfreqs2
 inquire(file='tcoord.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading tcoord from tcoord.txt'
  tcoord=0d0
  open(unit=iunit, form='formatted', file='tcoord.txt', status='old')
  read(iunit,*) ntimes_
  do itau=1,ntimes
   read(iunit,*) tcoord(itau)
  enddo
  close(iunit)
 endif
 inquire(file='tweight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading tweight from tweight.txt'
  tweight=0d0
  open(unit=iunit, form='formatted', file='tweight.txt', status='old')
  read(iunit,*) ntimes_
  do itau=1,ntimes
   read(iunit,*) tweight(itau)
  enddo
  close(iunit)
 endif
 inquire(file='wcoord.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading wcoord2 from wcoord.txt'
  wcoord2=0d0
  open(unit=iunit, form='formatted', file='wcoord.txt', status='old')
  read(iunit,*) nfreqs2_
  do jfreq=1,nfreqs2
   read(iunit,*) wcoord2(jfreq)
  enddo
  close(iunit)
 endif
 inquire(file='wweight.txt', exist=file_exists)
 if(file_exists) then
  write(*,*) 'Reading wweight2 from wweight.txt'
  wweight2=0d0
  open(unit=iunit, form='formatted', file='wweight.txt', status='old')
  read(iunit,*) nfreqs2_
  do jfreq=1,nfreqs2
   read(iunit,*) wweight2(jfreq)
  enddo
  close(iunit)
 endif
 write(*,*)
!-----------!
! scGW loop !
!-----------!

 iter=0
 do
  iter=iter+1
  ! For iter=1 we build G_ao_itau as the RHF one
  if(iter==1) then
   G_ao_itau=czero
   do itau=1,ntimes
    call G0itau_ao_RHF(nBas,nOrb,nO, tcoord(itau),G_plus_itau ,cHF,eHF)
    call G0itau_ao_RHF(nBas,nOrb,nO,-tcoord(itau),G_minus_itau,cHF,eHF)
    G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:)
    G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:)
   enddo
  endif

  ! Build using the time grid Xo(i tau) = -2i G(i tau) G(-i tau)
  !  then Fourier transform Xo(i tau) -> Xo(i w)
  Chi0_ao_iw(:,:,:)=czero
  do itau=1,ntimes
   ! Xo(i tau) = -2i G(i tau) G(-i tau)
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas                       
                                   ! r1   r2'                    r2   r1'
       product = G_ao_itau(2*itau-1,ibas,jbas)*G_ao_itau(2*itau,kbas,lbas)
       if(abs(product)<1e-12) product=czero
       Chi0_ao_itau(1+(lbas-1)+(ibas-1)*nBas,1+(kbas-1)+(jbas-1)*nBas) = product
      enddo
     enddo
    enddo
   enddo
   Chi0_ao_itau=-2d0*im*Chi0_ao_itau ! The 2 factor is added to account for both spin contributions [ i.e., (up,up,up,up) and (down,down,down,down) ]
   ! Xo(i tau) -> Xo(i w) [the factor fact cancells large itau and iw values that lead to large oscillations]
   do ifreq=1,nfreqs
    Chi0_ao_iw(ifreq,:,:) = Chi0_ao_iw(ifreq,:,:) - im*tweight(itau)*Chi0_ao_itau(:,:)*Exp(im*tcoord(itau)*wcoord(ifreq))
   enddo 
  enddo
   ! Complete the Xo(i tau) -> Xo(i w)
   Chi0_ao_iw(:,:,:) = 2d0*Real(Chi0_ao_iw(:,:,:))

  ! Build Wp_ao_iw(i w) and Sigma_c_ao(i w2)
  !    and compute Ec Galitskii-Migdal
  EcGM=0d0
  Sigma_c_ao=0d0
  do ifreq=1,nfreqs
   trace1=0d0; trace2=0d0;
   ! Xo(i w) -> Wp_ao_iw(i w)
   Wp_ao_iw(:,:)=-matmul(Real(Chi0_ao_iw(ifreq,:,:)),vMAT(:,:))  
   do ibas=1,nBas2
    trace1=trace1+Wp_ao_iw(ibas,ibas)
    Wp_ao_iw(ibas,ibas)=Wp_ao_iw(ibas,ibas)+1d0
   enddo
   call inverse_matrix(nBas2,Wp_ao_iw,Wp_ao_iw)
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),Real(Chi0_ao_iw(ifreq,:,:)))
   Wp_ao_iw(:,:)=matmul(Wp_ao_iw(:,:),vMAT(:,:))
   do ibas=1,nBas2
    trace2=trace2+Wp_ao_iw(ibas,ibas)
   enddo
   EcGM=EcGM-wweight(ifreq)*(trace2+trace1)/(2d0*pi) ! iw contribution to EcGM
   Wp_ao_iw(:,:)=matmul(vMAT(:,:),Wp_ao_iw(:,:)) ! Now Wp_ao_iw is on the iw grid
   ! Wp_ao_iw(i w) -> Sigma_c_ao(i w2) 
   do jfreq=1,nfreqs2
    ! Build G(i (w+w2)) and G(i (w-w2))
    weval(1)=wcoord2(jfreq)+wcoord(ifreq)
    weval(2)=wcoord2(jfreq)-wcoord(ifreq)
    ! Doing G(itau) -> G(iw) using Fourier transform
    call Gitau2Giw(nBas,ntimes,ntimes_twice,tweight,tcoord,weval,G_ao_itau,G_ao_1,G_ao_2)
!    if(iter==1) then ! For iter=1 we could use the analytic expresion of G_AO_RHF
!     weval_cpx=im*weval(1)
!     call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,G_ao_1)
!     weval_cpx=im*weval(2)
!     call G_AO_RHF(nBas,nOrb,nO,eta,cHF,eHF,weval_cpx,G_ao_2)
!    endif
    ! Build Sigma_c_ao(i w2)
    do ibas=1,nBas
     do jbas=1,nBas
      do kbas=1,nBas
       do lbas=1,nBas
        Sigma_c_ao(jfreq,ibas,jbas)=Sigma_c_ao(jfreq,ibas,jbas)-(G_ao_1(kbas,lbas)+G_ao_2(kbas,lbas)) &
                                   *Wp_ao_iw(1+(kbas-1)+(ibas-1)*nBas,1+(jbas-1)+(lbas-1)*nBas)       &
                                   *wweight(ifreq)/(2d0*pi)
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
  ! Converge with respect to the Fock operator (using only good P_ao matrices)
  iter_fock=0
  do
   ! Build F
   iter_fock=iter_fock+1
   F_ao=Hc
   Ehfl=0d0
   do ibas=1,nBas
    do jbas=1,nBas
     Ehfl=Ehfl+P_ao(ibas,jbas)*Hc(ibas,jbas)
     do kbas=1,nBas
      do lbas=1,nBas
       F_ao(ibas,jbas)=F_ao(ibas,jbas)+P_ao(kbas,lbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
                      -0.5d0*P_ao(kbas,lbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
       Ehfl=Ehfl+0.5d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(lbas-1)+(kbas-1)*nBas,1+(jbas-1)+(ibas-1)*nBas) &
           -0.25d0*P_ao(kbas,lbas)*P_ao(ibas,jbas)*vMAT(1+(jbas-1)+(kbas-1)*nBas,1+(lbas-1)+(ibas-1)*nBas)
      enddo
     enddo
    enddo
   enddo
   ! Build G(i w2) and n(r)
   P_ao_old=P_ao
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao,wcoord2,wweight2, &
                      G_ao_1,G_ao_iw2,P_ao,trace_1_rdm) 
   if(abs(trace_1_rdm-nElectrons)>thrs_N) &
    call fix_chem_pot_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao,wcoord2,wweight2, &
                          G_ao_1,G_ao_iw2,P_ao,trace_1_rdm)
   ! Check convergence of P_ao
   diff_Pao=0d0
   do ibas=1,nBas
    do jbas=1,nBas
     diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_old(ibas,jbas))
    enddo
   enddo
   if(diff_Pao<=thrs_Pao) exit

   if(iter_fock==maxSCF) exit

  enddo

  ! Check convergence
  diff_Pao=0d0
  do ibas=1,nBas
   do jbas=1,nBas
    diff_Pao=diff_Pao+abs(P_ao(ibas,jbas)-P_ao_iter(ibas,jbas))
   enddo
  enddo
  P_ao_iter=P_ao

  ! Print iter info
  P_mo=-matmul(matmul(cHFinv,P_ao),transpose(cHFinv)) ! Minus to order occ numbers
  call diagonalize_matrix(nOrb,P_mo,Occ)
  write(*,*)
  write(*,'(a,f15.8,a,i5,a,i5)') ' Trace scGW  ',trace_1_rdm,' after ',iter_fock,' Fock iterations at global iter ',iter
  write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
  write(*,'(a,f15.8)')        ' Chem. Pot.  ',chem_pot
  write(*,'(a,f15.8)')        ' EcGM        ',EcGM
  write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
  write(*,*)
  write(*,*) ' Occupation numbers'
  Occ=-Occ
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,Occ(ibas)
  enddo

  if(diff_Pao<=thrs_Pao) exit

  if(iter==maxSCF) exit

  ! Prepare the a and b coefs that fit G_pq(i w2_k) ~ f_pq(i w2_k) = a_pq / (i w2_k -b_pq )
  call fit_a_b_coefs(nBas,nfreqs2,wcoord2,G_ao_iw2,a_coef,b_coef) 

  ! Compute newG_pq(i w2_k) = G_pq(i w2_k) - f_pq(i -w2_k)
  do jfreq=1,nfreqs2
   do ibas=1,nBas
    do jbas=1,nBas
     G_ao_iw2(jfreq,ibas,jbas)=G_ao_iw2(jfreq,ibas,jbas)-a_coef(ibas,jbas)/(im*wcoord2(jfreq)-b_coef(ibas,jbas))
    enddo
   enddo
  enddo 

  ! Build G(i w2) -> G(i tau) [ i tau and -i tau ]
  !   newG_ao_iw2(nfreqs2) --> G_ao_itau(ntimes_twice) + f(i +-tau)
  G_ao_itau=czero
  do itau=1,ntimes
   call Giw2Gitau(nBas,nfreqs2,wweight2,wcoord2,tcoord(itau),G_ao_iw2,G_plus_itau,G_minus_itau)
   call build_tail_fpq_itau(nBas, tcoord(itau),a_coef,b_coef,f_pq_tau)
   G_ao_itau(2*itau-1,:,:)=G_plus_itau(:,:)  + f_pq_tau(:,:)
   call build_tail_fpq_itau(nBas,-tcoord(itau),a_coef,b_coef,f_pq_tau)
   G_ao_itau(2*itau  ,:,:)=G_minus_itau(:,:) + f_pq_tau(:,:)
  enddo

 enddo
 write(*,*)
 write(*,'(A50)') '---------------------------------------'
 write(*,'(A50)') '      scGW calculation completed       '
 write(*,'(A50)') '---------------------------------------'
 write(*,*)
 write(*,'(a,f15.8,a,i5,a)') ' Trace scGW  ',trace_1_rdm,' after ',iter,' global iterations '
 write(*,'(a,f15.8)')        ' Change of P ',diff_Pao
 write(*,'(a,f15.8)')        ' Chem. Pot.  ',chem_pot
 write(*,'(a,f15.8)')        ' Hcore+Hx    ',Ehfl
 write(*,'(a,f15.8)')        ' EcGM        ',EcGM
 write(*,'(a,f15.8)')        ' Eel         ',Ehfl+EcGM
 write(*,'(a,f15.8)')        ' Etot        ',Ehfl+EcGM+ENuc
 write(*,*)
 write(*,*) ' Final occupation numbers'
 do ibas=1,nOrb
  write(*,'(I7,F15.8)') ibas,Occ(ibas)
 enddo

 call wall_time(end_scGWitauiw)
 
 t_scGWitauiw = end_scGWitauiw - start_scGWitauiw
 write(*,'(A65,1X,F9.3,A8)') 'Total wall time for scGW = ',t_scGWitauiw,' seconds'
 write(*,*)

 ! Restore values and deallocate dyn arrays
 eHF(:) = eHF(:)+chem_pot_saved
 deallocate(Chi0_ao_iw,Wp_ao_iw)
 deallocate(tcoord,tweight) 
 deallocate(wcoord2,wweight2) 
 deallocate(G_ao_itau)
 deallocate(Sigma_c_ao,G_ao_iw2)
 deallocate(a_coef,b_coef,f_pq_tau) 
 deallocate(P_ao,P_ao_old,P_ao_iter,F_ao,P_mo,cHFinv,Occ) 
 deallocate(G_minus_itau,G_plus_itau) 
 deallocate(G_ao_1,G_ao_2) 
 deallocate(Chi0_ao_itau) 

end subroutine 

subroutine fix_chem_pot_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao,wcoord2,wweight2, &
                             G_ao,G_ao_iw2,P_ao,trace_1_rdm) 

! Fix the chemical potential for scGW 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs2
  double precision,intent(in)   :: nElectrons
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord2(nfreqs2)
  double precision,intent(in)   :: wweight2(nfreqs2)
  complex*16,intent(in)         :: Sigma_c_ao(nfreqs2,nBas,nBas)

! Local variables

  integer                       :: isteps
  double precision              :: thrs_closer
  double precision              :: delta_chem_pot
  double precision              :: chem_pot_change
  double precision              :: chem_pot_old
  double precision              :: grad_electrons
  double precision              :: trace_2up
  double precision              :: trace_up
  double precision              :: trace_down
  double precision              :: trace_2down
  double precision              :: trace_old

! Output variables

  double precision,intent(inout):: chem_pot
  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: P_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao_iw2(nfreqs2,nBas,nBas)

  !  Initialize 

  isteps = 0
  delta_chem_pot  = 2d-1
  thrs_closer     = 2d-1
  chem_pot_change = 0d0
  grad_electrons  = 1d0
  trace_1_rdm      = -1d0

  write(*,*)
  write(*,'(a)') ' Fixing the Tr[1D] at scGW '
  write(*,*)
  write(*,*)'------------------------------------------------------'
  write(*,'(1X,A1,1X,A15,1X,A1,1X,A15,1X,A1A15,2X,A1)') &
          '|','Tr[1D]','|','Chem. Pot.','|','Grad N','|'
  write(*,*)'------------------------------------------------------'

  ! First approach close the value with an error lower than 1

  trace_old = 1d2
  do while( abs(trace_old-nElectrons) > thrs_closer .and. isteps <= 100 )
   isteps = isteps + 1
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_old) 
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_down) 
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_up) 
   if( abs(trace_up-nElectrons) > abs(trace_old-nElectrons) .and. abs(trace_down-nElectrons) > abs(trace_old-nElectrons) ) then
     write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
     '|',trace_old,'|',chem_pot,'|',grad_electrons,'|'
     delta_chem_pot = 0.75d0*delta_chem_pot
     thrs_closer = 0.5d0*thrs_closer
     write(*,*) "| contracting ...                                     |"
     if(delta_chem_pot<1d-2) exit
   else
     if( abs(trace_up-nElectrons) < abs(trace_old-nElectrons) ) then
      chem_pot=chem_pot+delta_chem_pot
      write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
      '|',trace_up,'|',chem_pot,'|',grad_electrons,'|'
     else
      if( abs(trace_down-nElectrons) < abs(trace_old-nElectrons) ) then
       chem_pot=chem_pot-delta_chem_pot
       write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
       '|',trace_down,'|',chem_pot,'|',grad_electrons,'|'
      endif
     endif
   endif
  enddo

  ! Do  final search

  write(*,*)'------------------------------------------------------'
  isteps = 0
  delta_chem_pot  = 1.0d-3
  do while( abs(trace_1_rdm-nElectrons) > thrs_N .and. isteps <= 100 )
   isteps = isteps + 1
   chem_pot = chem_pot + chem_pot_change
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_1_rdm)
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot+2d0*delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_2up)
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot+delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_up)
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot-delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_down)
   call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot-2d0*delta_chem_pot,S,F_ao,Sigma_c_ao, &
                      wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_2down)
!   grad_electrons = (trace_up-trace_down)/(2d0*delta_chem_pot)
   grad_electrons = (-trace_2up+8d0*trace_up-8d0*trace_down+trace_2down)/(12d0*delta_chem_pot)
   chem_pot_change = -(trace_1_rdm-nElectrons)/(grad_electrons+1d-10)
   ! Maximum change is bounded within +/- 0.10
   chem_pot_change = max( min( chem_pot_change , 0.1d0 / real(isteps) ), -0.1d0 / real(isteps) )
   write(*,'(1X,A1,F16.10,1X,A1,F16.10,1X,A1F16.10,1X,A1)') &
   '|',trace_1_rdm,'|',chem_pot,'|',grad_electrons,'|'
  enddo
  write(*,*)'------------------------------------------------------'
  write(*,*)
  call get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao, &
                     wcoord2,wweight2,G_ao,G_ao_iw2,P_ao,trace_old) 

end subroutine

subroutine get_1rdm_scGW(nBas,nfreqs2,nElectrons,thrs_N,chem_pot,S,F_ao,Sigma_c_ao,wcoord2,wweight2, &
                         G_ao,G_ao_iw2,P_ao,trace_1_rdm) 

! Compute the scGW 1RDM

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs2
  double precision,intent(in)   :: chem_pot
  double precision,intent(in)   :: nElectrons
  double precision,intent(in)   :: thrs_N
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: F_ao(nBas,nBas)
  double precision,intent(in)   :: wcoord2(nfreqs2)
  double precision,intent(in)   :: wweight2(nfreqs2)
  complex*16,intent(in)         :: Sigma_c_ao(nfreqs2,nBas,nBas)

! Local variables

  integer                       :: jfreq
  integer                       :: ibas,jbas
  complex*16                    :: weval_cpx

! Output variables

  double precision,intent(out)  :: trace_1_rdm
  double precision,intent(out)  :: P_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao(nBas,nBas)
  complex*16,intent(out)        :: G_ao_iw2(nfreqs2,nBas,nBas)

  P_ao=0d0
  G_ao_iw2=0d0
  trace_1_rdm=0d0
  do jfreq=1,nfreqs2
   weval_cpx=im*wcoord2(jfreq)
   ! Setting G(w) = [ (w+chem_pot)S + F +Sigma_c ]^-1
   G_ao(:,:)= (weval_cpx + chem_pot)*S(:,:) - F_ao(:,:) - Sigma_c_ao(jfreq,:,:) ! Gnew(iw2)^-1
   call complex_inverse_matrix(nBas,G_ao,G_ao)                                  ! Gnew(iw2)
   G_ao_iw2(jfreq,:,:)=G_ao(:,:)
   P_ao(:,:) = P_ao(:,:) + wweight2(jfreq)*real(G_ao(:,:)+conjg(G_ao(:,:)))
  enddo
  P_ao(:,:) = P_ao(:,:)/pi + S(:,:)
  do ibas=1,nBas
   do jbas=1,nBas
    trace_1_rdm=trace_1_rdm+P_ao(ibas,jbas)*S(ibas,jbas)
   enddo
  enddo

end subroutine

subroutine fit_a_b_coefs(nBas,nfreqs2,wcoord2,G_ao_iw2,a_coef,b_coef) 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nfreqs2

  double precision,intent(in)   :: wcoord2(nfreqs2)

  complex*16,intent(in)         :: G_ao_iw2(nfreqs2,nBas,nBas)

! Local variables

  logical                       :: line_bracket
  logical                       :: line_stage1
  integer                       :: ibas,jbas
  integer                       :: ifreq
  integer                       :: icall
  integer                       :: iflag
  integer                       :: lbfgs_status
  integer                       :: ndim
  integer                       :: history_record
  integer                       :: iter 
  integer                       :: line_info
  integer                       :: line_infoc
  integer                       :: line_nfev
  integer                       :: nwork
  double precision              :: norm_Gpq
  double precision              :: gtol
  double precision              :: line_stp
  double precision              :: line_stpmin
  double precision              :: line_stpmax
  double precision              :: line_dginit
  double precision              :: line_finit
  double precision              :: line_stx
  double precision              :: line_fx
  double precision              :: line_dgx
  double precision              :: line_sty
  double precision              :: line_fy
  double precision              :: line_dgy
  double precision              :: line_stmin
  double precision              :: line_stmax
  double precision              :: a,b,wk
  double precision              :: error
  double precision              :: ReG,ImG
  double precision, allocatable :: diag(:)
  double precision, allocatable :: work(:)
  double precision              :: ab_val(2)
  double precision              :: ab_grad(2)

! Output variables

  double precision,intent(out)  :: a_coef(nBas,nBas)
  double precision,intent(out)  :: b_coef(nBas,nBas)

 ! Initialization

 lbfgs_status = 0
 ndim   = 2
 history_record = 5
 nwork = ndim * ( 2 * history_record + 1 ) + 2 * history_record
 gtol = 0.90d0
 line_stpmin = 1.0d-20
 line_stpmax = 1.0d20
 line_stp    = 1.0d0

 allocate(work(nwork),diag(ndim))
 a_coef=0d0; b_coef=0d0;

 ! Use L-BFGS to optimize
 do ibas=1,nBas
  do jbas=1,nBas

   ! Initialize
   iter=0
   icall=0
   iflag=0
   diag(:) = 1d0
   work(:) = 0d0
   ab_val(:) = 0d0
   norm_Gpq = 0d0
   ! Optimize
   do
    error=0d0
    ab_grad=0d0
    do ifreq=nfreqs2-19,nfreqs2 ! Use the last 20 points to fit the tail
     a=ab_val(1)
     b=ab_val(2)
     wk=wcoord2(ifreq)
     ReG= Real( G_ao_iw2(ifreq,ibas,jbas) )
     if(abs(ReG)<1d-12) ReG=0d0
     ImG=Aimag( G_ao_iw2(ifreq,ibas,jbas) )
     if(abs(ImG)<1d-12) ImG=0d0
     if(icall==0) norm_Gpq=norm_Gpq+ReG**2d0+ImG**2d0
     error=error+(  -a*b/(b*b+wk*wk) - ReG )**2d0
     error=error+( -a*wk/(b*b+wk*wk) - ImG )**2d0
     ab_grad(1)=ab_grad(1)+( -a*b/(b*b+wk*wk)  -  ReG ) * ( -b/(b*b+wk*wk)  )
     ab_grad(1)=ab_grad(1)+( -a*wk/(b*b+wk*wk) -  ImG ) * ( -wk/(b*b+wk*wk) )
     ab_grad(2)=ab_grad(2)+( -a*b/(b*b+wk*wk)  -  ReG ) * ( a*(b-wk)*(b+wk)/(b*b+wk*wk)**2d0 )
     ab_grad(2)=ab_grad(2)+( -a*wk/(b*b+wk*wk) -  ImG ) * ( 2d0*a*b*wk/(b*b+wk*wk)**2d0 )
    enddo
    ab_grad(:)=2d0*ab_grad(:)
    if(norm_Gpq < 1d-12) exit
    if(sum(abs(ab_grad(:))) < 1d-8) exit

    call lbfgs(ndim, history_record, ab_val, error, ab_grad, diag, work, lbfgs_status, &
         gtol, line_stpmin, line_stpmax, line_stp, iter, line_info, line_nfev,         &
         line_dginit, line_finit,line_stx,  line_fx,  line_dgx,                        &
         line_sty,  line_fy,  line_dgy, line_stmin,  line_stmax,                       &
         line_bracket, line_stage1, line_infoc)

    iflag = lbfgs_status

    if(iflag<=0) exit
    icall=icall+1
    !  We allow at most 2000 evaluations
    if(icall==2000) exit

   enddo
   if(norm_Gpq < 1d-12) then
    ab_val(1) = 0d0
    ab_val(2) = 1d0
   endif
   ! Save the result
   a_coef(ibas,jbas)=ab_val(1)
   b_coef(ibas,jbas)=ab_val(2)
  enddo
 enddo
  
 ! Deallocate and exit
 deallocate(work,diag)

end subroutine
   
subroutine build_tail_fpq_itau(nBas,tcoord,a,b,f_pq_tau)

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: tcoord
  double precision,intent(in)   :: a(nBas,nBas)
  double precision,intent(in)   :: b(nBas,nBas)

! Local variables 

  integer                       :: ibas,jbas
  double precision              :: factor
  double precision,external     :: Heaviside_step

! Output variables

  complex*16,intent(out)        :: f_pq_tau(nBas,nBas)
  
  f_pq_tau=czero
  do ibas=1,nBas 
   do jbas=1,nBas
    factor=(Heaviside_step(-b(ibas,jbas))*Heaviside_step(tcoord)-Heaviside_step(b(ibas,jbas))*Heaviside_step(-tcoord))
    f_pq_tau(ibas,jbas)=im*a(ibas,jbas)*Exp(b(ibas,jbas)*tcoord)*factor
   enddo
  enddo

end subroutine
