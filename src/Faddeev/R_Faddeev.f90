subroutine R_Faddeev(working_dir,backward,dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,eHF,ERI,dipole, &
                     maxSCF_CC,thresh_CC,max_diis_CC,TDA_W,lin_GW,eta_GW,reg_GW,reg,eta,eta_screen)
  use R_FRPA_operator, only: run_faddeev_pairs, ccd_sources
  implicit none
  include 'parameters.h'
  character(*),intent(in) :: working_dir
  logical,intent(in) :: backward,dotest,TDA_W,lin_GW,reg_GW,reg
  integer,intent(in) :: nBas,nOrb,nC,nO,nV,nR,nS,maxSCF_CC,max_diis_CC
  double precision,intent(in) :: ENuc,ERHF,eHF(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),dipole(nOrb,nOrb,ncart)
  double precision,intent(in) :: thresh_CC,eta_GW,eta,eta_screen
  integer :: reference,pairs,coupling,nroots
  logical :: static3,external_gw,iterative
  double precision :: qp(nOrb),ePA(nOrb),ZPA(nOrb)
  double precision,allocatable :: t(:,:,:,:),Teh(:,:),Thh(:,:),Tee(:,:),ccd_ph(:,:),ccd_pp(:,:)

  call read_faddeev_options(working_dir,reference,pairs,coupling,static3,external_gw,iterative,nroots)
  if(iterative.and.nroots>nO-nC) error stop 'Faddeev IP root count exceeds active occupied space'
  if(pairs==2.and.coupling==4) error stop 'Dressed sources for CCD pair interactions are not defined'
  if(coupling==5.or.pairs==2) then
    allocate(t(nO-nC,nO-nC,nV-nR,nV-nR))
    call R_CCD_amplitudes(dotest,maxSCF_CC,thresh_CC,max_diis_CC,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eHF,t)
  end if
  if(pairs==1) then
    if(.not.allocated(t)) allocate(t(1,1,1,1),source=0d0)
    call R_GW_Faddeev(dotest,TDA_W,lin_GW,eta_GW,reg_GW,backward,eta,eta_screen,reg, &
                       nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole,eHF,qp,ePA,ZPA, &
                       coupling,iterative,nroots,static3,external_gw,t)
    return
  end if

  qp=eHF
  if(reference==1.or.external_gw) then
    call RG0W0(dotest,.false.,.false.,.false.,.false.,.false.,TDA_W,.false., &
                .false.,.false.,.false.,.true.,.true.,lin_GW,eta_GW,reg_GW,.false., &
                nOrb,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole,eHF,qp)
  end if
  if(coupling==4) then
    call R_bare_dressed_sources(nOrb,nC,nO,nR,nS,qp,ERI,reg,eta,Teh,Thh,Tee)
  end if
  if(coupling==5.or.pairs==2) then
    call ccd_sources(t,Teh,Thh)
    Tee=Thh
    if(pairs==2) then
      if(backward) then
        write(*,*) 'Faddeev: CCD pair interactions with projected RPA metric'
      else
        write(*,*) 'Faddeev: static CCD pair interactions'
      end if
      ccd_ph=Teh; ccd_pp=Thh
    end if
    deallocate(t)
  end if
  if(reference==0.and.external_gw) then
    call run_faddeev_pairs(backward,nOrb,nC,nO,nR,eHF,eHF,ERI,coupling,static3,.false.,iterative,nroots, &
                          reg,eta,Teh=Teh,Thh=Thh,Tee=Tee,ePA=ePA,ZPA=ZPA,external_energies=qp,ccd_ph=ccd_ph,ccd_pp=ccd_pp)
  else
    call run_faddeev_pairs(backward,nOrb,nC,nO,nR,eHF,qp,ERI,coupling,static3,external_gw,iterative,nroots, &
                          reg,eta,Teh=Teh,Thh=Thh,Tee=Tee,ePA=ePA,ZPA=ZPA,ccd_ph=ccd_ph,ccd_pp=ccd_pp)
  end if
contains
  subroutine R_bare_dressed_sources(nOrb,nC,nO,nR,nS,e,ERI,reg,eta,Teh,Thh,Tee)
    integer,intent(in) :: nOrb,nC,nO,nR,nS
    double precision,intent(in) :: e(nOrb),ERI(nOrb,nOrb,nOrb,nOrb),eta
    logical,intent(in) :: reg
    double precision,allocatable,intent(out) :: Teh(:,:),Thh(:,:),Tee(:,:)
    type sources
      double precision,allocatable :: eh(:,:),hh(:,:),ee(:,:)
    end type
    type(sources) :: source(2)
    double precision,allocatable :: Aph(:,:),Bph(:,:),Cee(:,:),Dhh(:,:),Bpp(:,:)
    integer :: nocc,nvir,nh,ne,spin,ia,i,a
    double precision :: shift
    logical :: shifted
    nocc=nO-nC; nvir=nOrb-nO-nR
    inquire(file='use_chem_pot_in_ppRPA',exist=shifted)
    shift=0d0
    if(shifted) shift=e(nO)+e(nO+1)
    allocate(Aph(nS,nS),Bph(nS,nS))
    do spin=1,2
      nh=nocc*(nocc+3-2*spin)/2
      ne=nvir*(nvir+3-2*spin)/2
      allocate(Cee(ne,ne),Dhh(nh,nh),Bpp(ne,nh))
      allocate(source(spin)%eh(nS,nS),source(spin)%hh(nh,ne),source(spin)%ee(nh,ne))
      call phRLR_A(spin,.false.,nOrb,nC,nO,nOrb-nO,nR,nS,1d0,e,ERI,Aph)
      call phRLR_B(spin,.false.,nOrb,nC,nO,nOrb-nO,nR,nS,1d0,ERI,Bph)
      ia=0
      do i=nC+1,nO
        do a=nO+1,nOrb-nR
          ia=ia+1
          Aph(ia,ia)=Aph(ia,ia)-e(a)+e(i)
        end do
      end do
      call ppRLR_C(spin,nOrb,nC,nO,nOrb-nO,nR,ne,1d0,e,ERI,Cee)
      call ppRLR_D(spin,nOrb,nC,nO,nOrb-nO,nR,nh,1d0,e,ERI,Dhh)
      call ppRLR_B(spin,nOrb,nC,nO,nOrb-nO,nR,nh,ne,1d0,ERI,Bpp)
      call TDA_Faddeev_strip_pair(nOrb,nC+1,nO,spin==1,.true.,nh,e,shift,Dhh)
      call TDA_Faddeev_strip_pair(nOrb,nO+1,nOrb-nR,spin==1,.false.,ne,e,shift,Cee)
      call TDA_Faddeev_sources(reg,eta,nOrb,nC,nO,nR,nS,nh,ne,spin==1,e, &
                              Aph,Dhh,Cee,Bph,Bpp,source(spin)%eh,source(spin)%hh,source(spin)%ee)
      deallocate(Cee,Dhh,Bpp)
    end do
    allocate(Teh(4*nS,4*nS),Thh(nocc*(2*nocc-1),nvir*(2*nvir-1)),Tee(nocc*(2*nocc-1),nvir*(2*nvir-1)))
    call R_TDA_Faddeev_lift_sources(nOrb,nC,nO,nR,nS,nocc*(nocc+1)/2,nocc*(nocc-1)/2,nvir*(nvir+1)/2,nvir*(nvir-1)/2, &
         source(1)%eh,source(2)%eh,source(1)%hh,source(2)%hh,source(1)%ee,source(2)%ee,Teh,Thh,Tee)
  end subroutine
end subroutine
