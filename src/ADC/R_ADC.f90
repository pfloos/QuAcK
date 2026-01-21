subroutine R_ADC(dotest,                                               & 
                 do_IPEA_ADC2,do_IPEA_ADC3,                            &
                 do_SOSEX,do_2SOSEX,do_G3W2,                           &
                 do_ADC_GW,do_ADC_2SOSEX,                              &
                 do_ADC3_G3W2,do_ADC3x_G3W2,do_ADC4_G3W2,              &
                 TDA_W,TDA,singlet,triplet,linearize,eta,doSRG,        &
                 do_dyson,diag_approx,sig_inf,                         & 
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

  logical,intent(in)            :: do_SOSEX
  logical,intent(in)            :: do_2SOSEX
  logical,intent(in)            :: do_G3W2

  logical,intent(in)            :: do_ADC_GW
  logical,intent(in)            :: do_ADC_2SOSEX
  logical,intent(in)            :: do_ADC3_G3W2
  logical,intent(in)            :: do_ADC3x_G3W2
  logical,intent(in)            :: do_ADC4_G3W2

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG
  
  logical,intent(in)            :: do_dyson
  logical,intent(in)            :: diag_approx     
  logical,intent(in)            :: sig_inf

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
  double precision,parameter    :: flow = 1d6

  logical                       :: do_hierarchy_GW = .true.
  logical                       :: do_1h1p,do_1h,do_diag
  logical                       :: do_full_freq,do_half_half, do_pure_stat

! Output variables
  
  ! None

  do_IPEA = do_IPEA_ADC2 .or. do_IPEA_ADC3 .or.       &  
            do_SOSEX .or. do_2SOSEX .or. do_G3W2 .or. &
            do_ADC_GW .or. do_ADC_2SOSEX .or. do_ADC3_G3W2 .or. do_ADC3x_G3W2 .or. do_ADC4_G3W2

  do_EE   = .false.

!=========================================!
! Charged excitation branch of ADC module !
!=========================================!

  if(do_IPEA) then 

  !----------------------------------!
  ! Perform IP/EA-ADC(2) calculation !
  !----------------------------------!

    if(do_IPEA_ADC2) then 
      
      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          call R_IPEA_ADC2_diag(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          call R_IPEA_ADC2(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
          call R_IP_ADC2_diag(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          print*, 'Full version of IP-ADC(2) not yet implemented'
        end if
   
      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for IP/EA-ADC(2) = ',t_ADC,' seconds'
      write(*,*)
 
    end if

  !----------------------------------!
  ! Perform IP/EA-ADC(3) calculation !
  !----------------------------------!

    if(do_IPEA_ADC3) then

      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          print*, 'Diagonal version of IPEA-ADC(3) not yet implemented'
!         call R_IPEA_ADC3(dotest,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          print*, 'Full version of IPEA-ADC(3) not yet implemented'
        end if

      else

        if(diag_approx) then
          print*, 'Diagonal version of IP-ADC(3) not yet implemented'
        else
          print*, 'Full version of IP-ADC(3) not yet implemented'
        end if

      end if
      call wall_time(end_ADC)

      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for IP/EA-ADC(3) = ',t_ADC,' seconds'
      write(*,*)

    end if

  !----------------------------!
  ! Perform SOSEX calculation !
  !----------------------------!

    if(do_SOSEX) then 
      
      call wall_time(start_ADC)
      call R_SOSEX(dotest,TDA_W,singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for SOSEX = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !----------------------------!
  ! Perform 2SOSEX calculation !
  !----------------------------!

    if(do_2SOSEX) then 
      
      call wall_time(start_ADC)
      call R_2SOSEX(dotest,TDA_W,singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for 2SOSEX = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !--------------------------!
  ! Perform G3W2 calculation !
  !--------------------------!

    if(do_G3W2) then 
      
      call wall_time(start_ADC)
      call R_G3W2(dotest,TDA_W,singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int_MO,eHF)
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for G3W2 = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !----------------------------!
  ! Perform ADC-GW calculation !
  !----------------------------!

    do_1h1p = .true.
    do_1h   = .true.
    do_diag = .true.

    do_full_freq = .true.
    do_half_half = .true.
    do_pure_stat = .true.

    if(do_hierarchy_GW) then 
      

      call wall_time(start_ADC)
      if(do_1h1p) then
          if(do_full_freq) call R_ADC_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_half_half) call R_half_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_pure_stat) call R_static_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      end if

      if(do_1h) then
          if(do_full_freq) call R_IP_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_half_half) call R_IP_half_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_pure_stat) call R_IP_static_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      end if

      if(do_diag) then
          if(do_full_freq) call R_ADC_GW_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_half_half) call R_diag_half_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          if(do_pure_stat) call R_diag_static_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC-GW = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !----------------------------!
  ! Perform ADC-GW calculation !
  !----------------------------!

    if(do_ADC_GW) then 
      
      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          call R_ADC_GW_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
!         call R_ADC_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
!         call R_IP_ADC_GW_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
!         call R_IP_ADC_GW(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      end if
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
      if(do_dyson) then

        if(diag_approx) then
          call R_ADC_2SOSEX_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          call R_ADC_2SOSEX(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
          print*, 'Diagonal version of IP-ADC-2SOSEX not yet implemented'
        else
          call R_IP_ADC_2SOSEX(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC-2SOSEX = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !---------------------------------!
  ! Perform ADC(3)-G3W2 calculation !
  !---------------------------------!

    if(do_ADC3_G3W2) then 
      
      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          call R_ADC3_G3W2_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          call R_ADC3_G3W2(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
          print*, 'Diagonal version of IP-ADC(3)-G3W2 not yet implemented'
        else
          print*, 'Full version of IP-ADC(3)-G3W2 not yet implemented'
        end if

      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC(3)-G3W2 = ',t_ADC,' seconds'
      write(*,*)
 
    end if
  
  !----------------------------------!
  ! Perform ADC(3x)-G3W2 calculation !
  !----------------------------------!

    if(do_ADC3x_G3W2) then 
      
      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          call R_ADC3x_G3W2_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          call R_ADC3x_G3W2(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
          print*, 'Diagonal version of IP-ADC(3x)-G3W2 not yet implemented'
        else
          print*, 'Full version of IP-ADC(3x)-G3W2 not yet implemented'
        end if

      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC(3x)-G3W2 = ',t_ADC,' seconds'
      write(*,*)
 
    end if

  !---------------------------------!
  ! Perform ADC(4)-G3W2 calculation !
  !---------------------------------!

    if(do_ADC4_G3W2) then 
      
      call wall_time(start_ADC)
      if(do_dyson) then

        if(diag_approx) then
          call R_ADC4_G3W2_diag(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
          ! call R_ADC4_G3W2_diag_fullmat(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        else
          call R_ADC4_G3W2(dotest,sig_inf,TDA_W,flow,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF)
        end if

      else

        if(diag_approx) then
          print*, 'Diagonal version of IP-ADC(4)-G3W2 not yet implemented'
        else
          print*, 'Full version of IP-ADC(4)-G3W2 not yet implemented'
        end if

      end if
      call wall_time(end_ADC)
    
      t_ADC = end_ADC - start_ADC
      write(*,'(A65,1X,F9.3,A8)') 'Total wall time for ADC(4)-G3W2 = ',t_ADC,' seconds'
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
