subroutine dfRG0W0mat(nOrb,nO,eta,shift,eHF,vMAT,nfreqs,ntimes,wcoord,wweight)

! Restricted branch of G0W0 matrix form

   implicit none
   include 'parameters.h'

   integer,intent(in)           :: nOrb
   integer,intent(in)           :: nO
   integer,intent(in)           :: nfreqs
   integer,intent(in)           :: ntimes

   double precision,intent(in)  :: eta
   double precision,intent(in)  :: shift
   double precision,intent(in)  :: eHF(nOrb)
   double precision,intent(in)  :: wcoord(nfreqs)
   double precision,intent(in)  :: wweight(nfreqs)
   double precision,intent(in)  :: vMAT(nOrb*nOrb,nOrb*nOrb)

   integer                      :: iorb,iter,istate,niter
   integer                      :: nE_eval_global
                                   
   double precision             :: chem_pot
   double precision,allocatable :: eQP_state(:)
   double precision,allocatable :: Hmat(:,:)
   double precision,allocatable :: Eigvec(:,:)
   double precision,allocatable :: Sigc_mo(:,:)
   double precision,allocatable :: Sigc_mo_tmp(:,:,:)
                                   
   complex*16,allocatable       :: E_eval_global_cpx(:)
   complex*16,allocatable       :: Sigc_mo_he_cpx(:,:,:)

   nE_eval_global=2
   niter=10

   allocate(eQP_state(nOrb))
   allocate(Eigvec(nOrb,nOrb))
   allocate(Hmat(nOrb,nOrb))
   allocate(Sigc_mo(nOrb,nOrb))
   allocate(Sigc_mo_he_cpx(nE_eval_global,nOrb,nOrb))
   allocate(E_eval_global_cpx(nE_eval_global))

   chem_pot=0.5d0*(eHF(nO)+eHF(nO+1))

   write(*,*)
   write(*,*)'**************************************'
   write(*,*)'*  G0W0 beyond the linearized approx *'
   write(*,*)'*     and using imag. freqs.         *'
   write(*,*)'*    [ Fock + Sigma_c(w) ] C = w C   *'
   write(*,*)'**************************************'
   write(*,*)

   do istate=1,nOrb
    write(*,*)
    write(*,*) '******'
    write(*,*) 'State',istate
    write(*,*)
    eQP_state(:)=eHF(:)-chem_pot
    eQP_state(:)=eQP_state(istate)
    do iter=1,niter
     write(*,*) 'Iter ',iter
     ! Set eQP to evaluate Sigma_c
     E_eval_global_cpx(1)=eQP_state(istate)-shift
     E_eval_global_cpx(2)=eQP_state(istate)+shift
     call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eHF,nfreqs,ntimes,wweight,wcoord,vMAT,&
                             Sigc_mo_he_cpx)
     ! Sigma_c_he
     Sigc_mo(:,:) = 0.5d0*(Real(Sigc_mo_he_cpx(1,:,:))+Real(Sigc_mo_he_cpx(2,:,:)))
     write(*,*) 'Sigma_c'
     do iorb=1,nOrb
      write(*,'(*(f10.5))') Sigc_mo(iorb,:)
     enddo
     ! New eQP value
     Hmat=0d0
     Hmat(1:nOrb           ,1:nOrb           )=Sigc_mo(1:nOrb,1:nOrb)
     do iorb=1,nOrb
      Hmat(iorb,iorb)=Hmat(iorb,iorb)+eHF(iorb)-chem_pot
     enddo
     Eigvec=Hmat
     call diagonalize_matrix(nOrb,Eigvec,eQP_state)
     eQP_state(:)=eQP_state(istate)
     write(*,*) 'New QP energy'
     write(*,'(*(f10.5))') eQP_state(istate)
     write(*,*)
    enddo
   
    ! Print results
    write(*,'(A,f15.8)') ' QP energy (aligned incl. the chem pot) ',eQP_state(istate)+chem_pot
    ! Set eQP to evaluate Sigma_c
    E_eval_global_cpx(1)=eQP_state(istate)-shift
    E_eval_global_cpx(2)=eQP_state(istate)+shift
    call build_Sigmac_w_RHF(nOrb,nO,nE_eval_global,eta,0,E_eval_global_cpx,eHF,nfreqs,ntimes,wweight,wcoord,vMAT,&
                            Sigc_mo_he_cpx)
    ! Sigma_c
    Sigc_mo(:,:) = 0.5d0*(Real(Sigc_mo_he_cpx(1,:,:))+Real(Sigc_mo_he_cpx(2,:,:)))
    write(*,*) 'Sigma_c'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Sigc_mo(iorb,:)
    enddo
    ! New eQP value
    Hmat=0d0
    Hmat(1:nOrb           ,1:nOrb           )=Sigc_mo(1:nOrb,1:nOrb)
    do iorb=1,nOrb
     Hmat(iorb,iorb)=Hmat(iorb,iorb)+eHF(iorb)-chem_pot
    enddo
    write(*,*) 'H = F + Sigma_c(e_QP)'
    do iorb=1,nOrb
     write(*,'(*(f10.5))') Hmat(iorb,:)
    enddo
    Hmat=matmul(Hmat,Eigvec)
    Hmat=Hmat-eQP_state(istate)*Eigvec
    write(*,*) 'H C = e_QP C ?'
    write(*,'(*(f10.5))') Hmat(:,istate)
   enddo

   deallocate(Sigc_mo)
   deallocate(Sigc_mo_he_cpx)
   deallocate(eQP_state)
   deallocate(Eigvec)
   deallocate(Hmat)
   deallocate(E_eval_global_cpx)

end subroutine
