subroutine print_qsUGF2(nBas,nO,nSCF,Conv,thresh,eHF,eGF2,cGF2,PGF2,Ov,T,V,J,K, & 
                        ENuc,ET,EV,EJ,Ex,Ec,EqsGF2,SigC,Z,dipole)

! Print one-electron energies and other stuff for qsUGF2

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)                 :: nBas
  integer,intent(in)                 :: nO(nspin)
  integer,intent(in)                 :: nSCF
  double precision,intent(in)        :: ENuc
  double precision,intent(in)        :: ET(nspin)
  double precision,intent(in)        :: EV(nspin)
  double precision,intent(in)        :: EJ(nsp)
  double precision,intent(in)        :: Ex(nspin)
  double precision,intent(in)        :: Ec(nsp)
  double precision,intent(in)        :: EqsGF2
  double precision,intent(in)        :: Conv
  double precision,intent(in)        :: thresh
  double precision,intent(in)        :: eHF(nBas,nspin)
  double precision,intent(in)        :: eGF2(nBas,nspin)
  double precision,intent(in)        :: cGF2(nBas,nBas,nspin)
  double precision,intent(in)        :: PGF2(nBas,nBas,nspin)
  double precision,intent(in)        :: Ov(nBas,nBas)
  double precision,intent(in)        :: T(nBas,nBas)
  double precision,intent(in)        :: V(nBas,nBas)
  double precision,intent(in)        :: J(nBas,nBas,nspin)
  double precision,intent(in)        :: K(nBas,nBas,nspin)
  double precision,intent(in)        :: SigC(nBas,nBas,nspin)
  double precision,intent(in)        :: Z(nBas,nspin)
  double precision,intent(in)        :: dipole(ncart)

! Local variables

  integer                            :: p
  integer                            :: ispin,ixyz
  double precision                   :: HOMO(nspin)
  double precision                   :: LUMO(nspin)
  double precision                   :: Gap(nspin)
  double precision                   :: S_exact,S2_exact
  double precision                   :: S,S2
  double precision,external          :: trace_matrix

! HOMO and LUMO

  do ispin=1,nspin
    if(nO(ispin) > 0) then
      HOMO(ispin) = eGF2(nO(ispin),ispin)
      LUMO(ispin) = eGF2(nO(ispin)+1,ispin)
      Gap(ispin)  = LUMO(ispin) - HOMO(ispin)
    else
      HOMO(ispin) = 0d0
      LUMO(ispin) = eGF2(1,ispin)
      Gap(ispin)  = 0d0
    end if
  end do

  S2_exact = dble(nO(1) - nO(2))/2d0*(dble(nO(1) - nO(2))/2d0 + 1d0)
  S2 = S2_exact + nO(2) - sum(matmul(transpose(cGF2(:,1:nO(1),1)),matmul(Ov,cGF2(:,1:nO(2),2)))**2)

  S_exact = 0.5d0*dble(nO(1) - nO(2))
  S = -0.5d0 + 0.5d0*sqrt(1d0 + 4d0*S2)

! Dump results

  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  if(nSCF < 10) then
    write(*,'(1X,A21,I1,A2,A12)')'  Self-consistent qsG',nSCF,'F2',' calculation'
  else
    write(*,'(1X,A21,I2,A2,A12)')'  Self-consistent qsG',nSCF,'F2',' calculation'
  endif
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(A1,A3,A1,A30,A1,A30,A1,A30,A1,A30,A1)') &
            '|',' ','|','e_HF            ','|','Sig_c            ','|','Z            ','|','e_QP            ','|'
  write(*,'(A1,A3,A1,2A15,A1,2A15,A1,2A15,A1,2A15,A1)') &
            '|','#','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|','up     ','dw     ','|'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'

  do p=1,nBas
    write(*,'(A1,I3,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1,2F15.6,A1)') &
    '|',p,'|',eHF(p,1)*HaToeV,eHF(p,2)*HaToeV,'|',SigC(p,p,1)*HaToeV,SigC(p,p,2)*HaToeV,'|', &
              Z(p,1),Z(p,2),'|',eGF2(p,1)*HaToeV,eGF2(p,2)*HaToeV,'|'
  enddo

  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(2X,A10,I3)')   'Iteration ',nSCF
  write(*,'(2X,A14,F15.5)')'Convergence = ',Conv
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') 'qsUGF2 HOMO      energy:',maxval(HOMO(:))*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'qsUGF2 LUMO      energy:',minval(LUMO(:))*HaToeV,' eV'
  write(*,'(2X,A30,F15.6,A3)') 'qsUGF2 HOMO-LUMO gap   :',(minval(LUMO(:))-maxval(HOMO(:)))*HaToeV,' eV'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,'(2X,A30,F15.6,A3)') '    qsUGF2 total       energy:',ENuc + EqsGF2,' au'
  write(*,'(2X,A30,F15.6,A3)') '    qsUGF2 exchange    energy:',sum(Ex(:)),' au'
  write(*,'(2X,A30,F15.6,A3)') '    qsUGF2 correlation energy:',sum(Ec(:)),' au'
  write(*,*)'-------------------------------------------------------------------------------'&
              '------------------------------------------------'
  write(*,*)

! Dump results for final iteration

  if(Conv < thresh) then

    write(*,*)
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40)')              ' Summary              '
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron    energy: ',sum(ET(:))  + sum(EV(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron a  energy: ',ET(1) + EV(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' One-electron b  energy: ',ET(2) + EV(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic         energy: ',sum(ET(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic      a  energy: ',ET(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Kinetic      b  energy: ',ET(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Potential       energy: ',sum(EV(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Potential    a  energy: ',EV(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Potential    b  energy: ',EV(2),' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron    energy: ',sum(EJ(:)) + sum(Ex(:)) + sum(Ec(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron aa energy: ',EJ(1) + Ex(1) + Ec(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron ab energy: ',EJ(2) + Ec(2),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Two-electron bb energy: ',EJ(3) + Ex(2) + Ec(3),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Hartree         energy: ',sum(EJ(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      aa energy: ',EJ(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      ab energy: ',EJ(2),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Hartree      bb energy: ',EJ(3),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Exchange        energy: ',sum(Ex(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Exchange     a  energy: ',Ex(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Exchange     b  energy: ',Ex(2),' au'
    write(*,*)
    write(*,'(A40,1X,F16.10,A3)') ' Correlation     energy: ',sum(Ec(:)),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Correlation  aa energy: ',Ec(1),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Correlation  ab energy: ',Ec(2),' au'
    write(*,'(A40,1X,F16.10,A3)') ' Correlation  bb energy: ',Ec(3),' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,1X,F16.10,A3)') ' Electronic      energy: ',EqsGF2,' au'
    write(*,'(A40,1X,F16.10,A3)') ' Nuclear      repulsion: ',ENuc,' au'
    write(*,'(A40,1X,F16.10,A3)') ' qsUGF2          energy: ',ENuc + EqsGF2,' au'
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A40,F13.6)')        '  S (exact)          :',2d0*S_exact + 1d0
    write(*,'(A40,F13.6)')        '  S                  :',2d0*S       + 1d0
    write(*,'(A40,F13.6)')        ' <S**2> (exact)      :',S2_exact
    write(*,'(A40,F13.6)')        ' <S**2>              :',S2
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,'(A45)')              ' Dipole moment (Debye)    '
    write(*,'(19X,4A10)')         'X','Y','Z','Tot.'
    write(*,'(19X,4F10.6)')       (dipole(ixyz)*auToD,ixyz=1,ncart),norm2(dipole)*auToD
    write(*,'(A60)')              '-------------------------------------------------'
    write(*,*)

    ! Print orbitals

    write(*,'(A50)') '-----------------------------------------'
    write(*,'(A50)') 'qsUGF2 spin-up   orbital coefficients '
    write(*,'(A50)') '-----------------------------------------'
    call matout(nBas,nBas,cGF2(:,:,1))
    write(*,*)
    write(*,'(A50)') '-----------------------------------------'
    write(*,'(A50)') 'qsUGF2 spin-down orbital coefficients '
    write(*,'(A50)') '-----------------------------------------'
    call matout(nBas,nBas,cGF2(:,:,2))
    write(*,*)

  endif

end subroutine 
