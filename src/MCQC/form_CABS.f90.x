subroutine form_CABS(nBas_OBS,nBas_ABS,c_OBS,c_ABS,S_ABS)

! Perform configuration interaction single calculation`

  implicit none

! Input variables

  integer,intent(in)            :: nBas_OBS,nBas_ABS
  double precision,intent(in)   :: S_ABS(nBas,nBas),c_OBS(nBas_OBS,nBas_OBS)

! Local variables

  integer                       :: 
  double precision              :: thresh = 1d-07
  integer                       :: i,j,a,b

! Output variables

  double precision,intent(out)  :: c_ABS(nBas_ABS,nBas_ABS)

  allocate(c(nBas_ABS,nBAs_OBS))

  c = 0d0
  c(1:nBas_OBS,1:nBas_OBS) = c_OBS(1:nBas_OBS,1:nBAs_OBS)

  c_ABS = 0d0
  do i=1,nBas_ABS
    c_ABS(i,i) = 1d0
  enddo

  v_ABS = S_ABS

  call DiagMat(nBas_ABS,v_ABS,e_ABS)
  
  nLD = 0
  do i=1,nBas_ABS
    if(abs(e_ABS(i)) < thresh) nLD = nLD +1
  enddo
  write(*,*) 'Number of linear dependencies in ABS',nLD

  call DoSVD(nBas_ABS,S_ABS,u,v,w)

! do a SVD of S_ABS to get u, v and w

  X_ABS = 0d0
  do i=1,nBas_ABS
    do j=1,nBas_ABS
      do k=1,nBas_ABS
        X_ABS(i,k) = X_ABS(i,k) + v_ABS(i,j)*e_ABS(j)*v_ABS(k,j)
      enddo
    enddo
  enddo

  cp_ABS = matmul(X_ABS,c_ABS)

  S12 = matmul(transpose(c),matmul(S_ABS,cp_ABS))


end subroutine form_CABS
