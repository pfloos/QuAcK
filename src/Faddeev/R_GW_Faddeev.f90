subroutine R_GW_Faddeev(dotest,TDA_W,lin_GW,eta_GW,reg_GW,frpa,eta_1b,eta_2b,reg_PA, &
                        nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eGW,ePA,ZPA, &
                        faddeev_coupling,faddeev_davidson,n_IP_roots,faddeev_static3,gw_external,ccd_t)
  implicit none
  include 'parameters.h'
  logical,intent(in) :: dotest,TDA_W,lin_GW,reg_GW,frpa,reg_PA,faddeev_davidson,faddeev_static3,gw_external
  integer,intent(in) :: nOrb,nC,nO,nV,nR,nS,faddeev_coupling,n_IP_roots
  double precision,intent(in) :: eta_GW,eta_1b,eta_2b,ENuc,ERHF,ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in) :: dipole_int(nOrb,nOrb,ncart),eHF(nOrb)
  double precision,intent(in) :: ccd_t(merge(nO-nC,1,faddeev_coupling==5),merge(nO-nC,1,faddeev_coupling==5), &
                                      merge(nV-nR,1,faddeev_coupling==5),merge(nV-nR,1,faddeev_coupling==5))
  double precision,intent(out) :: eGW(nOrb),ePA(nOrb),ZPA(nOrb)
  type pair_interactions
    double precision,allocatable :: eh(:,:),hh(:,:),ee(:,:),Beh(:,:),Bpp(:,:)
  end type
  type(pair_interactions) :: pair(2)
  logical :: adc2,dressed,shifted_pp
  integer :: nOcc,nVir,nOOs,nVVs,nOOt,nVVt,spin,nOO,nVV,i,a,ia,info
  double precision :: pair_shift
  double precision,allocatable :: OmW(:),rhoW(:,:,:)
    double precision,allocatable :: Tehs(:,:),Teht(:,:),Thhs(:,:),Thht(:,:),Tees(:,:),Teet(:,:)
    double precision,allocatable :: Teh(:,:),Thh(:,:),Tee(:,:),ip_omega(:),ip_Z(:),ip_residual(:)
  nOcc=nO-nC; nVir=nOrb-nR-nO
  nOOs=nOcc*(nOcc+1)/2; nOOt=nOcc*(nOcc-1)/2
  nVVs=nVir*(nVir+1)/2; nVVt=nVir*(nVir-1)/2
  adc2=faddeev_coupling==2; dressed=faddeev_coupling==4
  if(reg_PA.and.eta_1b<=0d0) error stop 'GW/BSE requires a positive one-body SRG flow'
  allocate(OmW(nS),rhoW(nOrb,nOrb,nS))
  call R_GW_BSE_screening(dotest,TDA_W,lin_GW,eta_GW,reg_GW,reg_PA,eta_2b, &
                          nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eGW,OmW,rhoW)

  inquire(file='use_chem_pot_in_ppRPA',exist=shifted_pp)
  pair_shift=0d0
  if(shifted_pp) pair_shift=eGW(nO)+eGW(nO+1)
  do spin=1,2
    nOO=merge(nOOs,nOOt,spin==1); nVV=merge(nVVs,nVVt,spin==1)
    allocate(pair(spin)%eh(nS,nS),pair(spin)%Beh(nS,nS))
    allocate(pair(spin)%hh(nOO,nOO),pair(spin)%ee(nVV,nVV),pair(spin)%Bpp(nVV,nOO))
    call R_GW_BSE_eh_blocks(spin,.not.(frpa.or.dressed),nOrb,nC,nO,nV,nR,nS,eta_GW,eGW,ERI,OmW,rhoW, &
                            pair(spin)%eh,pair(spin)%Beh)
    ia=0
    do i=nC+1,nO
      do a=nO+1,nOrb-nR
        ia=ia+1
        pair(spin)%eh(ia,ia)=pair(spin)%eh(ia,ia)-eGW(a)+eGW(i)
      end do
    end do
    call R_GW_BSE_pp_blocks(spin,.not.(frpa.or.dressed),nOrb,nC,nO,nV,nR,nS,nOO,nVV,eta_GW,eGW,ERI,OmW,rhoW, &
                            pair(spin)%Bpp,pair(spin)%ee,pair(spin)%hh)
    call TDA_Faddeev_strip_pair(nOrb,nC+1,nO,spin==1,.true.,nOO,eGW,pair_shift,pair(spin)%hh)
    call TDA_Faddeev_strip_pair(nOrb,nO+1,nOrb-nR,spin==1,.false.,nVV,eGW,pair_shift,pair(spin)%ee)
  end do
  deallocate(OmW,rhoW)
    if(frpa) then
      write(*,*) 'GW-BSE FRPA: frozen GW internal energies'
    else if(gw_external) then
      write(*,*) 'GW-BSE FTDA: GW external block; frozen GW internal energies'
    else
      write(*,*) 'GW-BSE FTDA: HF external block; frozen GW internal energies'
    end if
    allocate(Teh(4*nS,4*nS),Thh(nOcc*(2*nOcc-1),nVir*(2*nVir-1)),Tee(nOcc*(2*nOcc-1),nVir*(2*nVir-1)))
    Teh=0d0
    Thh=0d0
    Tee=0d0
    if(dressed) then
      allocate(Tehs(nS,nS),Teht(nS,nS),Thhs(nOOs,nVVs),Thht(nOOt,nVVt),Tees(nOOs,nVVs),Teet(nOOt,nVVt))
      call TDA_Faddeev_sources(reg_PA,eta_1b,nOrb,nC,nO,nR,nS,nOOs,nVVs,.true.,eGW, &
           pair(1)%eh,pair(1)%hh,pair(1)%ee,pair(1)%Beh,pair(1)%Bpp,Tehs,Thhs,Tees)
      call TDA_Faddeev_sources(reg_PA,eta_1b,nOrb,nC,nO,nR,nS,nOOt,nVVt,.false.,eGW, &
           pair(2)%eh,pair(2)%hh,pair(2)%ee,pair(2)%Beh,pair(2)%Bpp,Teht,Thht,Teet)
      call R_TDA_Faddeev_lift_sources(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
           Tehs,Teht,Thhs,Thht,Tees,Teet,Teh,Thh,Tee)
    end if
    if(faddeev_coupling==5) then
      call R_Faddeev_ccd_sources(nOcc,nVir,ccd_t,Teh,Thh)
      Tee=Thh
      write(*,*) 'Faddeev: CCD-amplitude couplings on the HF reference'
    end if
    if(frpa) then
      call R_FRPA_BSE(nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt,eHF,eGW,ERI, &
                     pair(1)%eh,pair(2)%eh,pair(1)%hh,pair(2)%hh,pair(1)%ee,pair(2)%ee, &
                     pair(1)%Beh,pair(2)%Beh,pair(1)%Bpp,pair(2)%Bpp,Teh,Thh,Tee, &
                     faddeev_coupling,faddeev_static3,gw_external,faddeev_davidson,n_IP_roots,reg_PA,eta_1b,ePA,ZPA)
    else if(faddeev_davidson) then
      allocate(ip_omega(n_IP_roots),ip_Z(n_IP_roots),ip_residual(n_IP_roots))
      call R_TDA_Faddeev_IP_roots(n_IP_roots,adc2,reg_PA,eta_1b,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
           eHF,eGW,ERI,pair(1)%eh,pair(2)%eh,pair(1)%hh,pair(2)%hh,pair(1)%ee,pair(2)%ee, &
           ip_omega,ip_Z,ip_residual,info,dressed.or.faddeev_coupling==5,Teh,Thh,Tee,faddeev_static3,gw_external)
      if(info/=0) error stop 'GW-BSE FTDA Davidson failed to converge'
      ! Only the printed IP root table is a calculated spectrum in this mode.
      ePA=eHF
      ZPA=0d0
    else
      call R_TDA_Faddeev_self_energy(adc2,reg_PA,eta_1b,nOrb,nC,nO,nR,nS,nOOs,nOOt,nVVs,nVVt, &
           eHF,eGW,ERI,pair(1)%eh,pair(2)%eh,pair(1)%hh,pair(2)%hh,pair(1)%ee,pair(2)%ee, &
           ePA,ZPA,dressed.or.faddeev_coupling==5,Teh,Thh,Tee,faddeev_static3,gw_external)
    end if
end subroutine
