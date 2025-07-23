subroutine dfRG0W0Bmat(nBas,nO,nOrb,nOrb_twice,chem_pot,eta,shift,eQP_state,U_QP,vMAT,nfreqs,ntimes, &
                       wcoord,wweight,sigma,S,T,V,Hc,c,P,Panom,Delta,ERI)

! Restricted branch of G0W0 Bogoliubov matrix form

  implicit none
  include 'parameters.h'

  integer,intent(in)             :: nBas
  integer,intent(in)             :: nO
  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nOrb_twice
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: ntimes
                                 
  double precision,intent(inout) :: chem_pot
  double precision,intent(in)    :: sigma
  double precision,intent(in)    :: eta
  double precision,intent(in)    :: shift
  double precision,intent(in)    :: wcoord(nfreqs)
  double precision,intent(in)    :: wweight(nfreqs)
  double precision,intent(inout) :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)    :: vMAT(nOrb*nOrb,nOrb*nOrb)
  double precision,intent(in)    :: S(nBas,nBas)
  double precision,intent(in)    :: T(nBas,nBas)
  double precision,intent(in)    :: V(nBas,nBas)
  double precision,intent(in)    :: Hc(nBas,nBas) 
  double precision,intent(in)    :: c(nBas,nOrb)
  double precision,intent(in)    :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(inout) :: P(nBas,nBas)
  double precision,intent(inout) :: Panom(nBas,nBas)
  double precision,intent(inout) :: Delta(nBas,nBas)
                                 
  integer                        :: iorb,iter,istate,niter
  integer                        :: nE_eval_global
                                    
  double precision               :: trace_1rdm
  double precision               :: thrs_N
  double precision,allocatable   :: eQP_state_HFB(:)
  double precision,allocatable   :: R(:,:)
  double precision,allocatable   :: J(:,:)
  double precision,allocatable   :: K(:,:)
  double precision,allocatable   :: F(:,:)
  double precision,allocatable   :: Delta_MO(:,:)
  double precision,allocatable   :: Hmat(:,:)
  double precision,allocatable   :: Eigvec(:,:)
  double precision,allocatable   :: Sigc_mo_he(:,:)
  double precision,allocatable   :: Sigc_mo_hh(:,:)
  double precision,allocatable   :: Sigc_mo_eh(:,:)
  double precision,allocatable   :: Sigc_mo_ee(:,:)
  double precision,allocatable   :: Sigc_mo_tmp(:,:,:)
                                    
  complex*16,allocatable         :: E_eval_global_cpx(:)
  complex*16,allocatable         :: Sigc_mo_he_cpx(:,:,:)
  complex*16,allocatable         :: Sigc_mo_hh_cpx(:,:,:)
  complex*16,allocatable         :: Sigc_mo_eh_cpx(:,:,:)
  complex*16,allocatable         :: Sigc_mo_ee_cpx(:,:,:)

  double precision,intent(inout) :: eQP_state(nOrb_twice)

  nE_eval_global=2
  niter=10
  thrs_N=1d-8

  allocate(E_eval_global_cpx(nE_eval_global))
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(F(nOrb,nOrb))
  allocate(Delta_MO(nOrb,nOrb))
  allocate(eQP_state_HFB(nOrb_twice))
  allocate(Hmat(nOrb_twice,nOrb_twice))
  allocate(Eigvec(nOrb_twice,nOrb_twice))
  allocate(Sigc_mo_he(nOrb,nOrb))
  allocate(Sigc_mo_hh(nOrb,nOrb))
  allocate(Sigc_mo_eh(nOrb,nOrb))
  allocate(Sigc_mo_ee(nOrb,nOrb))
  allocate(Sigc_mo_he_cpx(nE_eval_global,nOrb,nOrb))
  allocate(Sigc_mo_hh_cpx(nE_eval_global,nOrb,nOrb))
  allocate(Sigc_mo_eh_cpx(nE_eval_global,nOrb,nOrb))
  allocate(Sigc_mo_ee_cpx(nE_eval_global,nOrb,nOrb))
  allocate(R(nOrb_twice,nOrb_twice))

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'*  G0W0 beyond the linearized approx *'
  write(*,*)'*     and using imag. freqs.         *'
  write(*,*)'*   [ H^HFB + Sigma_c(w) ] C = w C   *'
  write(*,*)'**************************************'
  write(*,*)

  eQP_state_HFB(:)=eQP_state(:)

  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)
  call anomalous_matrix_AO_basis(nBas,sigma,Panom,ERI,Delta)
  F=matmul(transpose(c),matmul(Hc+J+0.5d0*K-chem_pot*S,c))
  Delta_MO=matmul(transpose(c),matmul(Delta,c))
  ! F MO
  write(*,*)
  write(*,*) 'F'
  do iorb=1,nOrb
   write(*,'(*(f10.5))') F(iorb,:)
  enddo
  ! Delta MO
  write(*,*)
  write(*,*) 'Delta'
  do iorb=1,nOrb
   write(*,'(*(f10.5))') Delta(iorb,:)
  enddo
  write(*,*)

  do istate=1,nOrb_twice

   write(*,*)
   write(*,*) '******'
   write(*,*) 'State',istate
   write(*,*)
   eQP_state(:)=eQP_state_HFB(istate)
   do iter=1,niter
    write(*,*) 'Iter ',iter
    ! Set eQP to evaluate Sigma_c
    E_eval_global_cpx(1)=eQP_state(istate)-shift
    E_eval_global_cpx(2)=eQP_state(istate)+shift
    call build_Sigmac_w_RHFB(nOrb,nOrb_twice,nE_eval_global,eta,0,E_eval_global_cpx,eQP_state_HFB, &
                             nfreqs,ntimes,wweight,wcoord,vMAT,U_QP,Sigc_mo_he_cpx,Sigc_mo_hh_cpx, &
                             Sigc_mo_eh_cpx,Sigc_mo_ee_cpx,.true.,.true.)
    
    ! Sigma_c_he
    Sigc_mo_he(:,:) = 0.5d0*(Real(Sigc_mo_he_cpx(1,:,:))+Real(Sigc_mo_he_cpx(2,:,:)))
    write(*,*) 'Sigma_c,he'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Sigc_mo_he(iorb,:)
    enddo
    ! Sigma_c_hh
    Sigc_mo_hh(:,:) = 0.5d0*(Real(Sigc_mo_hh_cpx(1,:,:))+Real(Sigc_mo_hh_cpx(2,:,:)))
    write(*,*) 'Sigma_c,hh'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Sigc_mo_hh(iorb,:)
    enddo
    ! Sigma_c_eh
    Sigc_mo_eh(:,:) = 0.5d0*(Real(Sigc_mo_eh_cpx(1,:,:))+Real(Sigc_mo_eh_cpx(2,:,:)))
    write(*,*) 'Sigma_c,eh'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Sigc_mo_eh(iorb,:)
    enddo
    ! Sigma_c_ee
    Sigc_mo_ee(:,:) = 0.5d0*(Real(Sigc_mo_ee_cpx(1,:,:))+Real(Sigc_mo_ee_cpx(2,:,:)))
    write(*,*) 'Sigma_c,ee'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Sigc_mo_ee(iorb,:)
    enddo
    ! New eQP value
     ! Initialized H with Sigma_c^Gorkov. Note: We use Sigc_mo,hh/ee incl. a minus to recover the potential
    Hmat=0d0
    Hmat(1:nOrb           ,1:nOrb           )=Sigc_mo_he(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)=Sigc_mo_eh(1:nOrb,1:nOrb)
    Hmat(1:nOrb           ,nOrb+1:nOrb_twice)=-Sigc_mo_hh(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,1:nOrb           )=-Sigc_mo_ee(1:nOrb,1:nOrb)
     ! Add the Fock and Delta Contributions
    Hmat(1:nOrb           ,1:nOrb           )=Hmat(1:nOrb           ,1:nOrb           )+F(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)=Hmat(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)-F(1:nOrb,1:nOrb)
    Hmat(1:nOrb           ,nOrb+1:nOrb_twice)=Hmat(1:nOrb           ,nOrb+1:nOrb_twice)+Delta(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,1:nOrb           )=Hmat(nOrb+1:nOrb_twice,1:nOrb           )+Delta(1:nOrb,1:nOrb)
    write(*,*) 'H^HFB + Sigma_c(w) [ Note: using -Sigma_c,hh and -Sigma_c,ee ]'
    do iorb=1,nOrb_twice
     write(*,'(*(f10.5))') Hmat(iorb,:)
    enddo
     ! Diagonalize H^HFB + Sigma_c(w) 
    Eigvec=Hmat
    call diagonalize_matrix(nOrb_twice,Eigvec,eQP_state)
    ! Test R 
    R(:,:)=0d0
    do iorb=1,nOrb
     R(:,:)=R(:,:)+matmul(Eigvec(:,iorb:iorb),transpose(Eigvec(:,iorb:iorb))) 
    enddo
    trace_1rdm=0d0
    do iorb=1,nOrb
     trace_1rdm=trace_1rdm+R(iorb,iorb) 
    enddo
    write(*,'(a,f10.5)') ' Tr[^1D] =',2d0*trace_1rdm
    ! Adjusting the chemical potential is apparently not needed...
     !if( abs(trace_1rdm-nO) > thrs_N ) & 
     ! call fix_chem_pot(nO,nOrb,nOrb_twice,0,thrs_N,trace_1rdm,chem_pot,Hmat,Eigvec,R,&
     !                   eQP_state)
     ! U_QP=Eigvec 
    write(*,*) 'Eigenvalues'
    write(*,'(*(f10.5))') eQP_state(:)
    eQP_state(:)=eQP_state(istate)
    write(*,*) 'New QP energy'
    write(*,'(*(f10.5))') eQP_state(istate)
    write(*,*)
   enddo

   write(*,'(A,f15.8)') ' QP energy (aligned incl. the chem pot) ',eQP_state(istate)+chem_pot

  enddo
  
  deallocate(Sigc_mo_he)
  deallocate(Sigc_mo_hh)
  deallocate(Sigc_mo_eh)
  deallocate(Sigc_mo_ee)
  deallocate(Sigc_mo_he_cpx)
  deallocate(Sigc_mo_hh_cpx)
  deallocate(Sigc_mo_eh_cpx)
  deallocate(Sigc_mo_ee_cpx)
  deallocate(eQP_state_HFB)
  deallocate(Eigvec)
  deallocate(Hmat)
  deallocate(Delta_MO,F,J,K,R)
  deallocate(E_eval_global_cpx)


end subroutine
