subroutine R_ADC(dotest,                                               & 
                 do_IPEA_ADC2,do_IPEA_ADC3,                            &
                 do_ADC_GW,do_ADC_2SOSEX,do_ADC_G3W2,                  &
                 TDA_W,TDA,singlet,triplet,eta,doSRG,                  & 
                 nNuc,ZNuc,rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,         &
                 S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO, &
                 ERHF,PHF,FHF,cHF,eHF)

! Restricted ADC module

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: do_IPEA_ADC2
  logical,intent(in)            :: do_IPEA_ADC3
  logical,intent(in)            :: do_ADC_GW
  logical,intent(in)            :: do_ADC_2SOSEX
  logical,intent(in)            :: do_ADC_G3W2

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_MO(nOrb,nOrb,ncart)

  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: FHF(nBas,nBas)

! Local variables

  double precision              :: start_ADC,end_ADC,t_ADC
  logical                       :: do_IPEA,do_EE
  logical                       :: do_IP_ADC2 = .true.

! Output variables
  
  ! None

  do_IPEA = do_IPEA_ADC2 .or. do_IP_ADC2 .or. do_IPEA_ADC3 .or. & 
            do_ADC_GW .or. do_ADC_2SOSEX .or. do_ADC_G3W2
  do_EE   = .false.

!=========================================!
! Charged excitation branch of ADC module !
!=========================================!

  if(do_IPEA) then 

  !--------------------------------!
  ! Perform IP/EA-ADC2 calculation !
  !--------------------------------!

    if(do_IPEA_ADC2) then 
      
      call wall_time(start_ADC)
      call R_IPEA_ADC2(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for IP/EA-ADC(2) = ',t_ADC,' seconds'
      write(*,*)
 
    end if

  !-----------------------------------------!
  ! Perform (non-Dyson) IP-ADC2 calculation !
  !-----------------------------------------!
  
    if(do_IP_ADC2) then 
      
      call wall_time(start_ADC)
      call R_IP_ADC2(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for IP-ADC(2) = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !----------------------------!
  ! Perform ADC-GW calculation !
  !----------------------------!

    if(do_ADC_GW) then 
      
      call wall_time(start_ADC)
      call R_ADC_GW(dotest,TDA_W,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC-GW = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !--------------------------------!
  ! Perform ADC-2SOSEX calculation !
  !--------------------------------!

    if(do_ADC_2SOSEX) then 
      
      call wall_time(start_ADC)
      call R_ADC_2SOSEX(dotest,TDA_W,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC-2SOSEX = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !------------------------------!
  ! Perform ADC-G3W2 calculation !
  !------------------------------!

    if(do_ADC_G3W2) then 
      
      call wall_time(start_ADC)
      call R_ADC_G3W2(dotest,TDA_W,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC-G3W2 = ',t_ADC,' seconds'
      write(*,*)
 
    end if

  end if
  
!=========================================!
! Neutral excitation branch of ADC module !
!=========================================!

  if(do_EE) then 

  !--------------------------!
  ! Perform ADC2 calculation !
  !--------------------------!

!   if(do_EE_ADC2) then 
      
!     call wall_time(start_ADC)
!     call R_EE_ADC2()
!     call wall_time(end_ADC)
    
!     t_GW = end_GW - start_ADC2
!     write(*,'(A65,1X,F9.3,A8)') 'Total wall time for IP/EA-ADC2 = ',t_ADC,' seconds'
!     write(*,*)
 
!   end if

  end if

end subroutine
