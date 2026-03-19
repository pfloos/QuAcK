! See Phys. Rev. C 91, 064320 (2015)
!  Integrals transformations are given in the same order as in the paper Eqs. A1d to A1i

subroutine ERI_MO2QP_H40(nOrb2,ERI_MO_sw,Ua,Va,H_QP)
  implicit none
! Input variables
  integer,intent(in)            :: nOrb2
  double precision,intent(in)   :: Ua(nOrb2,nOrb2)
  double precision,intent(in)   :: Va(nOrb2,nOrb2)
  double precision,intent(in)   :: ERI_MO_sw(nOrb2,nOrb2,nOrb2,nOrb2)
! Local variables
  integer                       :: l1,l2,l3,l4
  integer                       :: k1,k2,k3,k4
  double precision,allocatable  :: TMP_MAT1(:,:,:,:)
  double precision,allocatable  :: TMP_MAT2(:,:,:,:)
! Output variables
  double precision,intent(inout) :: H_QP(nOrb2,nOrb2,nOrb2,nOrb2)
  allocate(TMP_MAT1(nOrb2,nOrb2,nOrb2,nOrb2))
  allocate(TMP_MAT2(nOrb2,nOrb2,nOrb2,nOrb2))

  TMP_MAT1=0d0
  do k1=1,nOrb2
   do l1=1,nOrb2
    do l2=1,nOrb2
     do l3=1,nOrb2
      do l4=1,nOrb2
        TMP_MAT1(k1,l2,l3,l4)=TMP_MAT1(k1,l2,l3,l4)+(ERI_MO_sw(l1,l2,l3,l4)-ERI_MO_sw(l1,l2,l4,l3))*Ua(l1,k1)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  TMP_MAT2=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do l2=1,nOrb2
     do l3=1,nOrb2
      do l4=1,nOrb2
        TMP_MAT2(k1,k2,l3,l4)=TMP_MAT2(k1,k2,l3,l4)+TMP_MAT1(k1,l2,l3,l4)*Ua(l2,k2)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 
  
  TMP_MAT1=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do k3=1,nOrb2
     do l3=1,nOrb2
      do l4=1,nOrb2
        TMP_MAT1(k1,k2,k3,l3)=TMP_MAT1(k1,k2,k3,l3)+TMP_MAT2(k1,k2,l3,l4)*Va(l4,k3)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  TMP_MAT2=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do k3=1,nOrb2
     do k4=1,nOrb2
      do l3=1,nOrb2
        H_QP(k1,k2,k3,k4)=H_QP(k1,k2,k3,k4)+TMP_MAT1(k1,k2,k3,l3)*Va(l3,k4)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  deallocate(TMP_MAT1,TMP_MAT2)
end subroutine

subroutine ERI_MO2QP_H40_2(nOrb2,ERI_MO_sw,Ua,Va,H_QP)
  implicit none
! Input variables
  integer,intent(in)            :: nOrb2
  double precision,intent(in)   :: Ua(nOrb2,nOrb2)
  double precision,intent(in)   :: Va(nOrb2,nOrb2)
  double precision,intent(in)   :: ERI_MO_sw(nOrb2,nOrb2,nOrb2,nOrb2)
! Local variables
  integer                       :: l1,l2,l3,l4
  integer                       :: k1,k2,k3,k4
  double precision,allocatable  :: TMP_MAT1(:,:,:,:)
  double precision,allocatable  :: TMP_MAT2(:,:,:,:)
! Output variables
  double precision,intent(inout) :: H_QP(nOrb2,nOrb2,nOrb2,nOrb2)
  allocate(TMP_MAT1(nOrb2,nOrb2,nOrb2,nOrb2))
  allocate(TMP_MAT2(nOrb2,nOrb2,nOrb2,nOrb2))

  TMP_MAT1=0d0
  do k1=1,nOrb2
   do l1=1,nOrb2
    do l2=1,nOrb2
     do l3=1,nOrb2
      do l4=1,nOrb2
        TMP_MAT1(k1,l2,l3,l4)=TMP_MAT1(k1,l2,l3,l4)-(ERI_MO_sw(l1,l2,l3,l4)-ERI_MO_sw(l1,l2,l4,l3))*Ua(l1,k1)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  TMP_MAT2=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do l2=1,nOrb2
     do l3=1,nOrb2
      do l4=1,nOrb2
        TMP_MAT2(k1,k2,l2,l3)=TMP_MAT2(k1,k2,l2,l3)+TMP_MAT1(k1,l2,l3,l4)*Va(l4,k2)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 
  
  TMP_MAT1=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do k3=1,nOrb2
     do l2=1,nOrb2
      do l3=1,nOrb2
        TMP_MAT1(k1,k2,k3,l3)=TMP_MAT1(k1,k2,k3,l3)+TMP_MAT2(k1,k2,l2,l3)*Ua(l2,k3)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  TMP_MAT2=0d0
  do k1=1,nOrb2
   do k2=1,nOrb2
    do k3=1,nOrb2
     do k4=1,nOrb2
      do l3=1,nOrb2
        H_QP(k1,k2,k3,k4)=H_QP(k1,k2,k3,k4)+TMP_MAT1(k1,k2,k3,l3)*Va(l3,k4)
      enddo 
     enddo 
    enddo 
   enddo 
  enddo 

  deallocate(TMP_MAT1,TMP_MAT2)
end subroutine


