subroutine s2_2rdm_HFB(nBas, nOrb, nOrb_twice, nO, Occ, sigma, cHFB, ERI)

! Print one-electron energies and other stuff

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas, nOrb, nOrb_twice
  integer,intent(in)                 :: nO
  double precision,intent(in)        :: sigma
  double precision,intent(in)        :: Occ(nOrb)
  double precision,intent(in)        :: cHFB(nBas,nOrb)
  double precision,intent(in)        :: ERI(nbas,nbas,nbas,nbas)

! Local variables

  integer                            :: iorb,jorb
  double precision                   :: Ne
  double precision                   :: Dijkl
  double precision                   :: Vee
  double precision                   :: Vee_aa
  double precision                   :: Vee_ab
  double precision                   :: trace_2rdm
  double precision                   :: trace_2rdm_aa
  double precision                   :: trace_2rdm_ab
  double precision                   :: s2_val
  double precision                   :: s2_val_aa
  double precision                   :: s2_val_ab

  double precision,allocatable       :: ERI_MO(:,:,:,:)

! Compute <S^2> using the NO 2-RDM of the HFB calculation
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,cHFB,ERI,ERI_MO)

  Vee=0d0
  Vee_aa=0d0
  Vee_ab=0d0
  s2_val_aa=0d0
  s2_val_ab=0d0
  trace_2rdm=0d0
  trace_2rdm_aa=0d0
  trace_2rdm_ab=0d0
  Ne=2d0*sum(Occ(:))
  s2_val=-Ne*(Ne-4d0)/4d0
  do iorb=1,nOrb
   do jorb=1,nOrb
    ! aa
     ! Hartree
    Dijkl=0.5d0*Occ(iorb)*Occ(jorb)    
    Vee_aa=Vee_aa+Dijkl*ERI_MO(iorb,jorb,iorb,jorb)
    trace_2rdm_aa=trace_2rdm_aa+Dijkl
    s2_val_aa=s2_val_aa+Dijkl 
     ! Exchange
    Dijkl=-0.5d0*Occ(iorb)*Occ(jorb)    
    Vee_aa=Vee_aa+Dijkl*ERI_MO(iorb,jorb,jorb,iorb)
    if(iorb==jorb) then
     trace_2rdm_aa=trace_2rdm_aa+Dijkl
     s2_val_aa=s2_val_aa+Dijkl 
    endif
    ! ab or ba
     ! Hartree
    Dijkl=0.5d0*Occ(iorb)*Occ(jorb)    
    Vee_ab=Vee_ab+Dijkl*ERI_MO(iorb,jorb,iorb,jorb)
    trace_2rdm_ab=trace_2rdm_ab+Dijkl
    if(iorb==jorb) then
     s2_val_ab=s2_val_ab+Dijkl
    endif
     ! Time-rev
    !Dijkl=-0.5d0*sqrt(Occ(iorb)*Occ(jorb)*(1d0-Occ(iorb))*(1d0-Occ(jorb))) ! CA NOFA v.2 (JKL)
    Dijkl=0.5d0*sigma*sqrt(Occ(iorb)*Occ(jorb)*(1d0-Occ(iorb))*(1d0-Occ(jorb)))
    Vee_ab=Vee_ab+Dijkl*ERI_MO(iorb,iorb,jorb,jorb)
    if(iorb==jorb) then
     trace_2rdm_ab=trace_2rdm_ab+Dijkl
     s2_val_ab=s2_val_ab+Dijkl
    endif
   enddo
  enddo
  trace_2rdm=2d0*(trace_2rdm_aa+trace_2rdm_ab)
  Vee=2d0*(Vee_aa+Vee_ab)
  s2_val=s2_val+2d0*(s2_val_aa-s2_val_ab)
  write(*,*) ' Quantities computed in MO basis for DEBUG'
  write(*,'(*(a,f17.8))') " Tr[ 2D   ] ",trace_2rdm
  write(*,'(*(a,f17.8))') " Tr[ 2Daa ] ",trace_2rdm_aa
  write(*,'(*(a,f17.8))') " Tr[ 2Dab ] ",trace_2rdm_ab
  write(*,'(*(a,f17.8))') " Vee        ",Vee
  write(*,'(*(a,f17.8))') " <S^2>      ",s2_val
  write(*,'(*(a,f17.8))') " <S^2>N     ",-Ne*(Ne-4d0)/4d0
  write(*,'(*(a,f17.8))') " <S^2>aa    ",2d0*s2_val_aa
  write(*,'(*(a,f17.8))') " <S^2>ab    ",-2d0*s2_val_ab
  deallocate(ERI_MO)

end subroutine 
