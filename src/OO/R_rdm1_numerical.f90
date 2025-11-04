subroutine R_rdm1_numerical(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,eHF,c,P,F)

       ! Calculate 1RDM numerically from RHF (or what ever you want)

  implicit none
  include 'parameters.h'

  ! Input
  logical,intent(in)            :: dotest
  logical,intent(in)            :: doaordm
  
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift
  
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  
  ! Local
  integer                       :: pind,qind,i
  double precision              :: factor
  double precision              :: Eplus
  double precision              :: Eminus
  double precision              :: delta = 1d-6
  
  double precision,allocatable  :: Hc_loc(:,:)
  double precision,allocatable  :: V_loc(:,:)
  double precision,allocatable  :: rdm1_MO(:,:)
  double precision,allocatable  :: rdm1(:,:)
  
  ! Output
    
  double precision,intent(out)  :: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)
                                                     
allocate(rdm1(nBas,nBas),rdm1_MO(nOrb,nOrb),Hc_loc(nBas,nBas),V_loc(nBas,nBas))

do pind=1,nBas
  do qind=pind,nBas
    Hc_loc(:,:) = Hc(:,:)
    V_loc(:,:)  = V(:,:)
    if(pind==qind) then                                    
      factor = 1d0                                   
      Hc_loc(pind,pind) = Hc(pind,pind) - delta                        
      V_loc(pind,pind)  =  V(pind,pind) - delta                        
    else                                             
      factor = 0.5d0                                 
      Hc_loc(pind,qind) = Hc(pind,qind) - delta                        
      Hc_loc(qind,pind) = Hc(qind,pind) - delta                        
      V_loc(pind,qind)  =  V(pind,qind) - delta                        
      V_loc(qind,pind)  =  V(qind,pind) - delta                        
    endif                                          
    call RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                 nBas,nOrb,nO,S,T,V_loc,Hc_loc,ERI,dipole_int,X,Eminus,eHF,c,P,F)
    
    ! Do eventually RPA                              
    Hc_loc(:,:) = Hc(:,:)
    if(pind==qind) then                                    
      factor = 1d0                                   
      Hc_loc(pind,pind) = Hc(pind,pind) + delta
      V_loc(pind,pind)  =  V(pind,pind) + delta                        
    else                                            
      factor = 0.5d0                                
      Hc_loc(pind,qind) = Hc(pind,qind) + delta                        
      Hc_loc(qind,pind) = Hc(qind,pind) + delta                        
      V_loc(pind,qind)  =  V(pind,qind) + delta                        
      V_loc(qind,pind)  =  V(qind,pind) + delta                        
    endif                                          
    call RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                 nBas,nOrb,nO,S,T,V_loc,Hc_loc,ERI,dipole_int,X,Eplus,eHF,c,P,F)
    ! Do eventually RPA                              
    
    ! Accounting for double occurence of offdiagonal terms
    rdm1(pind,qind) = factor*(Eplus - Eminus)/(2*delta)
    rdm1(qind,pind) = rdm1(pind,qind)
    write(*,*) Eplus, Eminus, Eplus - Eminus, delta
  enddo
enddo
call RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
              nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,Eminus,eHF,c,P,F)
write(*,*) "1RDM in AO"
call matout(nBas,nBas,rdm1)                                        
write(*,*) "1RDM in MO"
call AOtoMO(nBas,nOrb,c,rdm1,rdm1_MO)
call matout(nBas,nBas,rdm1_MO)                                        
deallocate(rdm1,rdm1_MO,Hc_loc,V_loc)

end subroutine
