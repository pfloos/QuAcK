
subroutine Sigma_c_GGF2B_brut(nBas2,nBas4,G_plus,G_minus,db_ERI_AO,Sigma_c_plus,Sigma_c_minus) 

! Restricted scGF2B all block M^8

  implicit none
  include 'parameters.h'

! Input variables
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nBas4

  double precision,intent(in)   :: db_ERI_AO(nBas2,nBas2,nBas2,nBas2)

  complex*16,intent(in)         :: G_plus(nBas4,nBas4)
  complex*16,intent(in)         :: G_minus(nBas4,nBas4)

! Local variables
  integer                       :: abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas

  complex*16                    :: Integ_val
  complex*16,allocatable        :: G_he_plus(:,:)
  complex*16,allocatable        :: G_hh_plus(:,:)
  complex*16,allocatable        :: G_ee_plus(:,:)
  complex*16,allocatable        :: G_eh_plus(:,:)
  complex*16,allocatable        :: G_he_minus(:,:)
  complex*16,allocatable        :: G_hh_minus(:,:)
  complex*16,allocatable        :: G_ee_minus(:,:)
  complex*16,allocatable        :: G_eh_minus(:,:)
  complex*16,allocatable        :: Sigma_he_plus(:,:)
  complex*16,allocatable        :: Sigma_hh_plus(:,:)
  complex*16,allocatable        :: Sigma_ee_plus(:,:)
  complex*16,allocatable        :: Sigma_eh_plus(:,:)
  complex*16,allocatable        :: Sigma_he_minus(:,:)
  complex*16,allocatable        :: Sigma_hh_minus(:,:)
  complex*16,allocatable        :: Sigma_ee_minus(:,:)
  complex*16,allocatable        :: Sigma_eh_minus(:,:)

! Output variables
  complex*16,intent(inout)      :: Sigma_c_plus(nBas4,nBas4)
  complex*16,intent(inout)      :: Sigma_c_minus(nBas4,nBas4)

  ! Allocate arrays

  allocate(G_he_plus(nBas2,nBas2)) 
  allocate(G_hh_plus(nBas2,nBas2))
  allocate(G_ee_plus(nBas2,nBas2))
  allocate(G_eh_plus(nBas2,nBas2))
  allocate(G_he_minus(nBas2,nBas2))
  allocate(G_hh_minus(nBas2,nBas2))
  allocate(G_ee_minus(nBas2,nBas2))
  allocate(G_eh_minus(nBas2,nBas2))
  allocate(Sigma_he_plus(nBas2,nBas2)) 
  allocate(Sigma_hh_plus(nBas2,nBas2))
  allocate(Sigma_ee_plus(nBas2,nBas2))
  allocate(Sigma_eh_plus(nBas2,nBas2))
  allocate(Sigma_he_minus(nBas2,nBas2))
  allocate(Sigma_hh_minus(nBas2,nBas2))
  allocate(Sigma_ee_minus(nBas2,nBas2))
  allocate(Sigma_eh_minus(nBas2,nBas2))

  ! Initialize 
  
  G_he_plus(1:nBas2,1:nBas2)=G_plus(      1:nBas2,      1:nBas2)
  G_hh_plus(1:nBas2,1:nBas2)=G_plus(      1:nBas2,nBas2+1:nBas4)    
  G_ee_plus(1:nBas2,1:nBas2)=G_plus(nBas2+1:nBas4,      1:nBas2)   
  G_eh_plus(1:nBas2,1:nBas2)=G_plus(nBas2+1:nBas4,nBas2+1:nBas4)  
  G_he_minus(1:nBas2,1:nBas2)=G_minus(      1:nBas2,      1:nBas2)
  G_hh_minus(1:nBas2,1:nBas2)=G_minus(      1:nBas2,nBas2+1:nBas4)    
  G_ee_minus(1:nBas2,1:nBas2)=G_minus(nBas2+1:nBas4,      1:nBas2)   
  G_eh_minus(1:nBas2,1:nBas2)=G_minus(nBas2+1:nBas4,nBas2+1:nBas4)  
  Sigma_he_plus=czero
  Sigma_hh_plus=czero
  Sigma_ee_plus=czero
  Sigma_eh_plus=czero
  Sigma_he_minus=czero
  Sigma_hh_minus=czero
  Sigma_ee_minus=czero
  Sigma_eh_minus=czero

  ! M^8 loop
  !write(*,*)
  !write(*,*) 'Computing Bog. Sigma_c (M^8)'
  !write(*,*)

  !$OMP PARALLEL DEFAULT(NONE)                                                    &
  !$OMP          PRIVATE(abas,bbas,cbas,dbas,gbas,ebas,fbas,hbas,Integ_val)       &
  !$OMP          SHARED(nBas2,nBas4,G_plus,G_minus,db_ERI_AO,                    &
  !$OMP &        G_he_plus,G_hh_plus,G_ee_plus,G_eh_plus,G_he_minus,G_hh_minus,   &
  !$OMP &        G_ee_minus,G_eh_minus,Sigma_he_plus,Sigma_hh_plus,Sigma_ee_plus, &
  !$OMP &        Sigma_eh_plus,Sigma_he_minus,Sigma_hh_minus,Sigma_ee_minus,Sigma_eh_minus)
  block ! Use block to define local arrays
  complex*16,allocatable        :: Sigma_he_plus_th(:,:)
  complex*16,allocatable        :: Sigma_hh_plus_th(:,:)
  complex*16,allocatable        :: Sigma_ee_plus_th(:,:)
  complex*16,allocatable        :: Sigma_eh_plus_th(:,:)
  complex*16,allocatable        :: Sigma_he_minus_th(:,:)
  complex*16,allocatable        :: Sigma_hh_minus_th(:,:)
  complex*16,allocatable        :: Sigma_ee_minus_th(:,:)
  complex*16,allocatable        :: Sigma_eh_minus_th(:,:)

  ! Allocate local arrays
  allocate(Sigma_he_plus_th(nBas2,nBas2)) 
  allocate(Sigma_hh_plus_th(nBas2,nBas2))
  allocate(Sigma_ee_plus_th(nBas2,nBas2))
  allocate(Sigma_eh_plus_th(nBas2,nBas2))
  allocate(Sigma_he_minus_th(nBas2,nBas2))
  allocate(Sigma_hh_minus_th(nBas2,nBas2))
  allocate(Sigma_ee_minus_th(nBas2,nBas2))
  allocate(Sigma_eh_minus_th(nBas2,nBas2))

  ! Initialize local arrays
  Sigma_he_plus_th=czero
  Sigma_hh_plus_th=czero
  Sigma_ee_plus_th=czero
  Sigma_eh_plus_th=czero
  Sigma_he_minus_th=czero
  Sigma_hh_minus_th=czero
  Sigma_ee_minus_th=czero
  Sigma_eh_minus_th=czero

  !$OMP DO COLLAPSE(1)
  do abas=1,nBas2
   do bbas=1,nBas2
    do cbas=1,nBas2
     do dbas=1,nBas2
      do ebas=1,nBas2
       do fbas=1,nBas2
        do gbas=1,nBas2
         do hbas=1,nBas2
           ! he 2'
           Integ_val=0.5d0*db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
           Sigma_he_plus_th(abas,bbas) = Sigma_he_plus_th(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_he_plus(fbas,gbas) *G_he_minus(hbas,ebas)            
           Sigma_he_minus_th(abas,bbas)=Sigma_he_minus_th(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_he_minus(fbas,gbas)*G_he_plus(hbas,ebas)            
           ! he 2''
           Integ_val=-db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(dbas,gbas,bbas,hbas)
           Sigma_he_plus_th(abas,bbas) = Sigma_he_plus_th(abas,bbas)+Integ_val*G_he_plus(cbas,dbas) *G_hh_plus(fbas,hbas) *G_ee_minus(gbas,ebas)            
           Sigma_he_minus_th(abas,bbas)=Sigma_he_minus_th(abas,bbas)+Integ_val*G_he_minus(cbas,dbas)*G_hh_minus(fbas,hbas)*G_ee_plus(gbas,ebas) 
           ! eh 2'
           Integ_val=0.5d0*db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(bbas,hbas,dbas,gbas)
           Sigma_eh_plus_th(abas,bbas) = Sigma_eh_plus_th(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_eh_plus(fbas,gbas) *G_eh_minus(hbas,ebas)            
           Sigma_eh_minus_th(abas,bbas)=Sigma_eh_minus_th(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_eh_minus(fbas,gbas)*G_eh_plus(hbas,ebas)            
           ! eh 2''
           Integ_val=-db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(bbas,hbas,dbas,gbas)
           Sigma_eh_plus_th(abas,bbas) = Sigma_eh_plus_th(abas,bbas)+Integ_val*G_eh_plus(cbas,dbas) *G_hh_plus(ebas,gbas) *G_ee_minus(hbas,fbas)            
           Sigma_eh_minus_th(abas,bbas)=Sigma_eh_minus_th(abas,bbas)+Integ_val*G_eh_minus(cbas,dbas)*G_hh_minus(ebas,gbas)*G_ee_plus(hbas,fbas) 
           ! hh 2'
           Integ_val=-db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
           Sigma_hh_plus_th(abas,bbas) = Sigma_hh_plus_th(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_he_plus(fbas,hbas) *G_he_minus(gbas,ebas)            
           Sigma_hh_minus_th(abas,bbas)=Sigma_hh_minus_th(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_he_minus(fbas,hbas)*G_he_plus(gbas,ebas) 
           ! hh 2''
           Integ_val=0.5d0*db_ERI_AO(abas,ebas,cbas,fbas)*db_ERI_AO(hbas,bbas,gbas,dbas)
           Sigma_hh_plus_th(abas,bbas) = Sigma_hh_plus_th(abas,bbas)+Integ_val*G_hh_plus(cbas,dbas) *G_hh_plus(fbas,gbas) *G_ee_minus(hbas,ebas)            
           Sigma_hh_minus_th(abas,bbas)=Sigma_hh_minus_th(abas,bbas)+Integ_val*G_hh_minus(cbas,dbas)*G_hh_minus(fbas,gbas)*G_ee_plus(hbas,ebas)            
           ! ee 2'
           Integ_val=-db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(gbas,dbas,hbas,bbas)
           Sigma_ee_plus_th(abas,bbas) = Sigma_ee_plus_th(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_he_plus(ebas,gbas) *G_he_minus(hbas,fbas)            
           Sigma_ee_minus_th(abas,bbas)=Sigma_ee_minus_th(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_he_minus(ebas,gbas)*G_he_plus(hbas,fbas) 
           ! ee 2''
           Integ_val=0.5*db_ERI_AO(cbas,fbas,abas,ebas)*db_ERI_AO(gbas,dbas,hbas,bbas)
           Sigma_ee_plus_th(abas,bbas) = Sigma_ee_plus_th(abas,bbas)+Integ_val*G_ee_plus(cbas,dbas) *G_hh_plus(ebas,hbas) *G_ee_minus(gbas,fbas)            
           Sigma_ee_minus_th(abas,bbas)=Sigma_ee_minus_th(abas,bbas)+Integ_val*G_ee_minus(cbas,dbas)*G_hh_minus(ebas,hbas)*G_ee_plus(gbas,fbas)            
         enddo
        enddo
       enddo
      enddo
     enddo
    enddo
   enddo
  enddo
  !$OMP END DO

  ! Add contributions (one thread at a time) 
  !$OMP CRITICAL

  Sigma_he_plus = Sigma_he_plus + Sigma_he_plus_th    
  Sigma_hh_plus = Sigma_hh_plus + Sigma_hh_plus_th
  Sigma_ee_plus = Sigma_ee_plus + Sigma_ee_plus_th
  Sigma_eh_plus = Sigma_eh_plus + Sigma_eh_plus_th
  Sigma_he_minus = Sigma_he_minus + Sigma_he_minus_th
  Sigma_hh_minus = Sigma_hh_minus + Sigma_hh_minus_th
  Sigma_ee_minus = Sigma_ee_minus + Sigma_ee_minus_th
  Sigma_eh_minus = Sigma_eh_minus + Sigma_eh_minus_th

  !$OMP END CRITICAL

  ! Wait for all to finish
  !$OMP BARRIER 

  ! Deallocate local arrays 
  deallocate(Sigma_he_plus_th) 
  deallocate(Sigma_hh_plus_th)
  deallocate(Sigma_ee_plus_th)
  deallocate(Sigma_eh_plus_th)
  deallocate(Sigma_he_minus_th)
  deallocate(Sigma_hh_minus_th)
  deallocate(Sigma_ee_minus_th)
  deallocate(Sigma_eh_minus_th)

  end block !
  !$OMP END PARALLEL

  ! Set Sigma_c Gorkov

  Sigma_c_plus(      1:nBas2,       1:nBas2)=Sigma_he_plus(1:nBas2,1:nBas2)
  Sigma_c_plus(      1:nBas2, nBas2+1:nBas4)=Sigma_hh_plus(1:nBas2,1:nBas2)    
  Sigma_c_plus(nBas2+1:nBas4,       1:nBas2)=Sigma_ee_plus(1:nBas2,1:nBas2)   
  Sigma_c_plus(nBas2+1:nBas4, nBas2+1:nBas4)=Sigma_eh_plus(1:nBas2,1:nBas2)  
  Sigma_c_minus(     1:nBas2,       1:nBas2)=Sigma_he_minus(1:nBas2,1:nBas2)
  Sigma_c_minus(     1:nBas2, nBas2+1:nBas4)=Sigma_hh_minus(1:nBas2,1:nBas2)    
  Sigma_c_minus(nBas2+1:nBas4,      1:nBas2)=Sigma_ee_minus(1:nBas2,1:nBas2)   
  Sigma_c_minus(nBas2+1:nBas4,nBas2+1:nBas4)=Sigma_eh_minus(1:nBas2,1:nBas2)  

  ! Deallocate arrays 

  deallocate(G_he_plus)
  deallocate(G_hh_plus)
  deallocate(G_ee_plus)
  deallocate(G_eh_plus)
  deallocate(G_he_minus)
  deallocate(G_hh_minus)
  deallocate(G_ee_minus)
  deallocate(G_eh_minus)
  deallocate(Sigma_he_plus) 
  deallocate(Sigma_hh_plus)
  deallocate(Sigma_ee_plus)
  deallocate(Sigma_eh_plus)
  deallocate(Sigma_he_minus)
  deallocate(Sigma_hh_minus)
  deallocate(Sigma_ee_minus)
  deallocate(Sigma_eh_minus)

end subroutine

