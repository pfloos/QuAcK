subroutine dfRG0W0Bmat(nOrb,nOrb_twice,eta,shift,eQP_state,U_QP,vMAT,nfreqs,ntimes,wcoord,wweight)

! Restricted branch of G0W0 Bogoliubov matrix form

  implicit none
  include 'parameters.h'

  integer,intent(in)             :: nOrb
  integer,intent(in)             :: nOrb_twice
  integer,intent(in)             :: nfreqs
  integer,intent(in)             :: ntimes
                                 
  double precision,intent(in)    :: eta
  double precision,intent(in)    :: shift
  double precision,intent(in)    :: wcoord(nfreqs)
  double precision,intent(in)    :: wweight(nfreqs)
  double precision,intent(in)    :: U_QP(nOrb_twice,nOrb_twice)
  double precision,intent(in)    :: vMAT(nOrb*nOrb,nOrb*nOrb)
                                 
  integer                        :: iorb,iter,istate,niter
  integer                        :: nE_eval_global
                                    
  double precision,allocatable   :: eQP_state_HFB(:)
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
  allocate(E_eval_global_cpx(nE_eval_global))

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'*  G0W0 beyond the linearized approx *'
  write(*,*)'*     and using imag. freqs.         *'
  write(*,*)'*   [ H^HFB + Sigma_c(w) ] C = w C   *'
  write(*,*)'**************************************'
  write(*,*)

  eQP_state_HFB(:)=eQP_state(:)

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
    Hmat=0d0
    Hmat(1:nOrb           ,1:nOrb           )=Sigc_mo_he(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,nOrb+1:nOrb_twice)=Sigc_mo_eh(1:nOrb,1:nOrb)
    Hmat(nOrb+1:nOrb_twice,1:nOrb           )=Sigc_mo_ee(1:nOrb,1:nOrb)
    Hmat(1:nOrb           ,nOrb+1:nOrb_twice)=Sigc_mo_hh(1:nOrb,1:nOrb)
    do iorb=1,nOrb_twice
     Hmat(iorb,iorb)=Hmat(iorb,iorb)+eQP_state_HFB(iorb)
    enddo
    Eigvec=Hmat
    call diagonalize_matrix(nOrb_twice,Eigvec,eQP_state)
    eQP_state(:)=eQP_state(istate)
    write(*,*) 'New QP energy'
    write(*,'(*(f10.5))') eQP_state(istate)
    write(*,*)
   enddo
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
  deallocate(E_eval_global_cpx)


end subroutine
