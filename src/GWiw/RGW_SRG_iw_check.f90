subroutine RGW_SRG_iw_check(nOrb,nC,nO,nV,nR,nS,eHF,ERI,nfreqs)

! Fixed-orbital check of the SRG static self-energy computed from imaginary-axis data only
! (RGW_SRG_iw_self_energy) against the pole-based reference of Marie and Loos (RGW_SRG_self_energy),
! both evaluated with the same RHF orbitals and orbital energies.

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS,nfreqs
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)

! Local variables

  integer,parameter             :: ns_flow = 4
  double precision,parameter    :: flows(ns_flow) = (/ 10d0, 100d0, 500d0, 1000d0 /)
  integer                       :: is,p,q
  double precision              :: EcRPA,EcGM,gap,Omin,Omax,emax,emin,errmax,errdiag
  double precision,allocatable  :: Aph(:,:),Bph(:,:),Om(:),XpY(:,:),XmY(:,:),rho(:,:,:)
  double precision,allocatable  :: SigRef(:,:),SigIw(:,:),Z(:)

  write(*,*)
  write(*,*) '**************************************************************************'
  write(*,*) '* SRG static self-energy: imaginary-axis version vs pole-based reference *'
  write(*,*) '**************************************************************************'
  write(*,*)

  allocate(Aph(nS,nS),Bph(nS,nS),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS))
  allocate(SigRef(nOrb,nOrb),SigIw(nOrb,nOrb),Z(nOrb))

! Reference screening: RPA excitation energies and spectral weights

  call phRLR_A(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR_B(1,.true.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
  call phRLR(.false.,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

  gap  = eHF(nO+1) - eHF(nO)
  emin = minval(eHF)
  emax = maxval(eHF)
  Omin = 0.5d0*gap
  Omax = 1.3d0*(emax-emin)
  write(*,'(A,2F12.5)') '   Lowest / highest RPA excitation energy       = ',Om(1),Om(nS)
  write(*,'(A,2F12.5)') '   Quadrature window [Omin,Omax] used (no poles) = ',Omin,Omax
  if(Om(1) < Omin .or. Om(nS) > Omax) write(*,*) '   WARNING: RPA poles fall outside the quadrature window'
  write(*,*)

  do is=1,ns_flow

    call RGW_SRG_self_energy(flows(is),nOrb,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigRef,Z)
    call RGW_SRG_iw_self_energy(flows(is),nOrb,nC,nO,nV,nR,eHF,ERI,nfreqs,SigIw)

    errmax  = maxval(abs(SigIw-SigRef))
    errdiag = 0d0
    do p=1,nOrb
      errdiag = max(errdiag,abs(SigIw(p,p)-SigRef(p,p)))
    end do

    write(*,'(A,F10.2)') '   SRG flow parameter s = ',flows(is)
    write(*,'(A)') '   orbital   e_HF          St_ref(p,p)       St_iw(p,p)        difference'
    do p=1,nOrb
      write(*,'(4X,I5,2X,F12.6,3X,F16.10,2X,F16.10,2X,ES12.3)') p,eHF(p),SigRef(p,p),SigIw(p,p),SigIw(p,p)-SigRef(p,p)
    end do
    write(*,'(A,ES12.3)') '   max |St_iw - St_ref| (diagonal)       = ',errdiag
    write(*,'(A,ES12.3)') '   max |St_iw - St_ref| (full matrix)    = ',errmax
    write(*,'(A,ES12.3)') '   max |St_ref| (full matrix)            = ',maxval(abs(SigRef))
    write(*,*)

  end do

  deallocate(Aph,Bph,Om,XpY,XmY,rho,SigRef,SigIw,Z)

end subroutine
