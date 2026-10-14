subroutine G_ADC_Parquet(nOrb,nC,nO,nV,nR,eHF,eQP,ERI,eh_Phi,pp_Phi,eADC,ZADC)

! Hermitian / PSD ADC-Parquet in the generalized spin-orbital formalism.
!
! Sigma_ADC-PA(w) = U^T [ w - (K + C_PA) ]^{-1} U
!
! This is a PSD ADC completion of the static parquet kernels, not an
! exact rewriting of the original parquet self-energy.
!
! At eh_Phi = pp_Phi = 0, this routine reduces exactly to the
! first-order U/K/C content of G_IPEA_ADC3.

  implicit none
  include 'parameters.h'

  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQP(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_Phi(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: pp_Phi(nOrb,nOrb,nOrb,nOrb)

  double precision,intent(out)  :: eADC(nOrb)
  double precision,intent(out)  :: ZADC(nOrb)

  integer                       :: nOcc,nVir,n2h1p,n2p1h,nH,nAct
  integer                       :: p,i,j,k,l,a,b,c,d,s
  integer                       :: ija,klc,iab,kcd,idx,it
  double precision              :: val,cutoff,asym,sumZ
  double precision              :: tol_deg,wgrp,wbest
  double precision              :: start_timing,end_timing,timing
  double precision,allocatable  :: H(:,:),ePA(:),Zstate(:)

  nOcc = nO - nC
  nVir = nOrb - nR - nO
  nAct = nOrb - nR - nC

  n2h1p = nOcc*(nOcc-1)*nVir/2
  n2p1h = nOcc*nVir*(nVir-1)/2
  nH    = nOrb + n2h1p + n2p1h

  cutoff  = 1d-2
  tol_deg = 1d-8

  allocate(H(nH,nH),ePA(nH),Zstate(nH))

  H(:,:)    = 0d0
  ePA(:)    = 0d0
  Zstate(:) = 0d0
  eADC(:)   = eHF(:)
  ZADC(:)   = 1d0

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Generalized Hermitian ADC-Parquet    *'
  write(*,*)'****************************************'
  write(*,*)
  write(*,'(1X,A20,I8)') '2h1p dimension     =',n2h1p
  write(*,'(1X,A20,I8)') '2p1h dimension     =',n2p1h
  write(*,'(1X,A20,I8)') 'Total dimension    =',nH
  write(*,*)

  call wall_time(start_timing)

! F block

  do p=1,nOrb
    H(p,p) = eHF(p)
  end do

! U^- and K^- : 2h1p

  ija = nOrb
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ija = ija + 1
        H(ija,ija) = eQP(i) + eQP(j) - eQP(a)

        do p=1,nOrb
          val = ERI(i,j,p,a) - ERI(i,j,a,p)
          H(p,ija) = val
          H(ija,p) = val
        end do

      end do
    end do
  end do

! C^- : exact index pattern of G_IPEA_ADC3 with
! bare pair interactions replaced by Gamma^eh / Gamma^pp

  ija = nOrb
  do i=nC+1,nO
    do j=i+1,nO
      do a=nO+1,nOrb-nR

        ija = ija + 1
        klc = nOrb

        do k=nC+1,nO
          do l=k+1,nO
            do c=nO+1,nOrb-nR

              klc = klc + 1
              val = 0d0

              if(a == c) val = val - Gamma_pp(i,j,k,l)
              if(i == k) val = val + Gamma_eh(c,j,a,l)
              if(j == l) val = val + Gamma_eh(c,i,a,k)
              if(i == l) val = val - Gamma_eh(c,j,a,k)
              if(j == k) val = val - Gamma_eh(c,i,a,l)

              H(ija,klc) = H(ija,klc) + val

            end do
          end do
        end do

      end do
    end do
  end do

! U^+ and K^+ : 2p1h

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        iab = iab + 1
        H(iab,iab) = eQP(a) + eQP(b) - eQP(i)

        do p=1,nOrb
          val = ERI(a,b,p,i) - ERI(a,b,i,p)
          H(p,iab) = val
          H(iab,p) = val
        end do

      end do
    end do
  end do

! C^+ : exact index pattern of G_IPEA_ADC3

  iab = nOrb + n2h1p
  do i=nC+1,nO
    do a=nO+1,nOrb-nR
      do b=a+1,nOrb-nR

        iab = iab + 1
        kcd = nOrb + n2h1p

        do k=nC+1,nO
          do c=nO+1,nOrb-nR
            do d=c+1,nOrb-nR

              kcd = kcd + 1
              val = 0d0

              if(i == k) val = val + Gamma_pp(a,b,c,d)
              if(a == c) val = val - Gamma_eh(k,b,i,d)
              if(b == d) val = val - Gamma_eh(k,a,i,c)
              if(a == d) val = val + Gamma_eh(k,b,i,c)
              if(b == c) val = val + Gamma_eh(k,a,i,d)

              H(iab,kcd) = H(iab,kcd) + val

            end do
          end do
        end do

      end do
    end do
  end do

! Check and enforce Hermiticity

  asym = maxval(abs(H(:,:) - transpose(H(:,:))))
  write(*,'(1X,A32,1X,ES12.4)') 'Max pre-symmetrization |H-H^T| =',asym

  H(:,:) = 0.5d0*(H(:,:) + transpose(H(:,:)))

  call wall_time(end_timing)
  timing = end_timing - start_timing

  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for ADC-PA matrix construction =',timing,' seconds'
  write(*,*)

! Diagonalize

  call wall_time(start_timing)
  call diagonalize_matrix(nH,H,ePA)
  call wall_time(end_timing)

  timing = end_timing - start_timing
  write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for ADC-PA diagonalization =',timing,' seconds'
  write(*,*)

! Spectral weights

  Zstate(:) = 0d0
  do s=1,nH
    do p=nC+1,nOrb-nR
      Zstate(s) = Zstate(s) + H(p,s)**2
    end do
  end do

  sumZ = sum(Zstate)
  write(*,'(1X,A30,F16.10)') 'Sum of active 1p weights      =',sumZ
  write(*,'(1X,A30,I16)')    'Expected active 1p dimension =',nAct
  write(*,*)

  do p=nC+1,nOrb-nR

    wbest = -1d0
    idx   = 1
    s     = 1

    do while(s <= nH)

      it = s
      do while(it < nH)
        if(abs(ePA(it+1) - ePA(s)) > tol_deg) exit
        it = it + 1
      end do

      wgrp = sum(H(p,s:it)**2)

      if(wgrp > wbest) then
        wbest = wgrp
        idx   = s
      end if

      s = it + 1

    end do

    eADC(p) = ePA(idx)
    ZADC(p) = wbest

  end do

! Print spectrum

  write(*,*)'-------------------------------------------------------'
  write(*,*)'         Generalized PSD ADC-Parquet spectrum           '
  write(*,*)'-------------------------------------------------------'
  write(*,'(1X,A6,2X,A16,2X,A16)') '#','energy (eV)','1p weight'

  do s=1,nH
    if(Zstate(s) > cutoff) then
      write(*,'(1X,I6,2X,F16.8,2X,F16.8)') s,ePA(s)*HaToeV,Zstate(s)
    end if
  end do

  write(*,*)'-------------------------------------------------------'
  write(*,*)

  deallocate(H,ePA,Zstate)

contains

  double precision function Gamma_eh(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_eh = ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)  &
             - eh_Phi(p1,p2,p4,p3)                   &
             + pp_Phi(p1,p2,p3,p4)

  end function Gamma_eh


  double precision function Gamma_pp(p1,p2,p3,p4)

    implicit none
    integer,intent(in) :: p1,p2,p3,p4

    Gamma_pp = ERI(p1,p2,p3,p4) - ERI(p1,p2,p4,p3)  &
             + eh_Phi(p1,p2,p3,p4)                   &
             - eh_Phi(p1,p2,p4,p3)

  end function Gamma_pp

end subroutine G_ADC_Parquet
