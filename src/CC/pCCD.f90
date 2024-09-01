
! ---

subroutine pCCD(dotest, maxIt, thresh, max_diis, nBas, nOrb, &
                nC, nO, nV, nR, Hc, ERI_AO, ENuc, ERHF, eHF, cHF)

! pair CCD module

  implicit none

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxIt
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas, nOrb, nC, nO, nV, nR
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: p,q,r,s,t,u,w
  integer                       :: pq,rs
  integer                       :: i,j,a,b

  integer                       :: nItAmp
  integer                       :: nItOrb
  double precision              :: CvgAmp
  double precision              :: CvgOrb
  double precision              :: ECC
  double precision              :: EcCC

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OV(:,:)

  double precision,allocatable  :: OOOO(:,:)
  double precision,allocatable  :: OOVV(:,:)
  double precision,allocatable  :: OVOV(:,:)
  double precision,allocatable  :: OVVO(:,:)
  double precision,allocatable  :: VVVV(:,:)

  double precision,allocatable  :: yO(:,:)
  double precision,allocatable  :: yV(:,:)

  double precision,allocatable  :: r2(:,:)
  double precision,allocatable  :: t2(:,:)
  double precision,allocatable  :: z2(:,:)

  double precision,allocatable  :: rdm1(:,:)
  double precision,allocatable  :: rdm2(:,:,:,:)
  double precision,allocatable  :: xOO(:,:)
  double precision,allocatable  :: xVV(:,:)
  double precision,allocatable  :: xOV(:,:)
  double precision              :: tr_1rdm
  double precision              :: tr_2rdm

  double precision              :: E1,E2
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: h(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)
  double precision,allocatable  :: grad(:)
  double precision,allocatable  :: tmp(:,:,:,:)
  double precision,allocatable  :: hess(:,:)
  double precision,allocatable  :: hessInv(:,:)
  double precision,allocatable  :: Kap(:,:)
  double precision,allocatable  :: ExpKap(:,:)

  integer                       :: O,V,N
  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: t2_diis(:,:)
  double precision,allocatable  :: z2_diis(:,:)
  double precision,external     :: trace_matrix
  double precision,external     :: Kronecker_delta
          
! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted pCCD Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Useful quantities

  O = nO - nC
  V = nV - nR
  N = O + V ! nOrb - nC - nR

  !------------------------------------!
  ! Star Loop for orbital optimization !
  !------------------------------------!

  allocate(ERI_MO(N,N,N,N))
  allocate(c(nBas,N), h(N,N))
  allocate(eO(O), eV(V), delta_OV(O,V))
  allocate(OOOO(O,O), OOVV(O,V), OVOV(O,V), OVVO(O,V), VVVV(V,V))

  do i = 1, N
    c(:,i) = cHF(:,nC+i)
  enddo

  CvgOrb = 1d0
  nItOrb = 0

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| Orbital Optimization for pCCD                     |'
  write(*,*)'----------------------------------------------------'

  do while(CvgOrb > thresh .and. nItOrb < 1)

    nItOrb = nItOrb + 1

    ! Transform integrals

    h = matmul(transpose(c), matmul(Hc, c))

    call AOtoMO_ERI_RHF(nBas, N, c(1,1), ERI_AO(1,1,1,1), ERI_MO(1,1,1,1))

    ! Form energy denominator

    eO(:) = eHF(nC+1:nO)
    eV(:) = eHF(nO+1:nOrb-nR)

    do i=1,O
      do a=1,V
        delta_OV(i,a) = eV(a) - eO(i)
      end do
    end do

    ! Create integral batches

    do i=1,O
      do j=1,O
        OOOO(i,j) = ERI_MO(i,i,j,j)
      end do
    end do
 
    do i=1,O
      do a=1,V
        OOVV(i,a) = ERI_MO(i,i,O+a,O+a)
        OVOV(i,a) = ERI_MO(i,O+a,i,O+a)
        OVVO(i,a) = ERI_MO(i,O+a,O+a,i)
      end do
    end do
 
    do a=1,V
      do b=1,V
        VVVV(a,b) = ERI_MO(O+a,O+a,O+b,O+b)
      end do
    end do

    !----------------------------!
    ! Star Loop for t amplitudes !
    !----------------------------!
 
    allocate(t2(O,V),r2(O,V),yO(O,O))
    allocate(err_diis(O*V,max_diis),t2_diis(O*V,max_diis))

    CvgAmp = 1d0
    nItAmp = 0
    ECC    = ERHF
    EcCC   = 0d0
 
    n_diis        = 0
    t2(:,:)       = 0d0
    t2_diis(:,:)  = 0d0
    err_diis(:,:) = 0d0
 
    write(*,*)
    write(*,*)'----------------------------------------------------'
    write(*,*)'| pCCD calculation: t amplitudes                   |'
    write(*,*)'----------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
              '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
    write(*,*)'----------------------------------------------------'
 
    do while(CvgAmp > thresh .and. nItAmp < maxIt)
 
    ! Increment 
 
      nItAmp = nItAmp + 1
 
    ! Form intermediate array
      
     yO(:,:) = matmul(t2,transpose(OOVV))
      
     ! Compute residual
 
      r2(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*t2(:,:) & 
              - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - OOVV(:,:)*t2(:,:))*t2(:,:)
 
      do i=1,O
        do a=1,V
 
          do j=1,O
            r2(i,a) = r2(i,a) - 2d0*OOVV(j,a)*t2(j,a)*t2(i,a) + OOOO(j,i)*t2(j,a) + yO(i,j)*t2(j,a) 
          end do 
 
          do b=1,V
            r2(i,a) = r2(i,a) - 2d0*OOVV(i,b)*t2(i,b)*t2(i,a) + VVVV(a,b)*t2(i,b)
          end do 
 
        end do
      end do
 
      ! Check convergence 
 
      CvgAmp = maxval(abs(r2(:,:)))
    
      ! Update amplitudes
 
      t2(:,:) = t2(:,:) - 0.5d0*r2(:,:)/delta_OV(:,:)
 
      ! Compute correlation energy
 
      EcCC = 0d0
      do i=1,O
        do a=1,V
          EcCC = EcCC + OOVV(i,a)*t2(i,a) 
        end do
      end do
 
      ! Dump results
 
      ECC = ERHF + EcCC
 
      ! DIIS extrapolation
 
      if(max_diis > 1) then
 
        n_diis = min(n_diis+1,max_diis)
        call DIIS_extrapolation(rcond,nO*nV,nO*nV,n_diis,err_diis,t2_diis,-0.5d0*r2/delta_OV,t2)
 
      end if
 
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
        '|',nItAmp,'|',ECC+ENuc,'|',EcCC,'|',CvgAmp,'|'
 
    end do
    write(*,*)'----------------------------------------------------'

    !---------------------------!
    ! End Loop for t amplitudes !
    !---------------------------!
 
    deallocate(r2,yO)
    deallocate(err_diis,t2_diis)

    ! Did it actually converge?
 
    if(nItAmp == maxIt) then
 
      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'! Convergence failed for t ampitudes !'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 
      stop
  
    end if
 
    !-----------------------------!
    ! Start Loop for z amplitudes !
    !-----------------------------!
 
    allocate(z2(O,V),r2(O,V),yO(O,O),yV(V,V))
    allocate(err_diis(O*V,max_diis),z2_diis(O*V,max_diis))

    CvgAmp = 1d0
    nItAmp = 0
 
    n_diis        = 0
    z2_diis(:,:)  = 0d0
    err_diis(:,:) = 0d0

    write(*,*)
    write(*,*)'----------------------------------------------------'
    write(*,*)'| pCCD calculation: z amplitudes                   |'
    write(*,*)'----------------------------------------------------'
    write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
              '|','#','|','E(pCCD)','|','Ec(pCCD)','|','Conv','|'
    write(*,*)'----------------------------------------------------'
 
    do while(CvgAmp > thresh .and. nItAmp < maxIt)
 
    ! Increment 
 
      nItAmp = nItAmp + 1
 
    ! Form intermediate array
    
     yO(:,:) = matmul(OOVV,transpose(t2))
     yV(:,:) = matmul(transpose(OOVV),t2)
 
     ! Compute residual
 
      r2(:,:) = OOVV(:,:) + 2d0*delta_OV(:,:)*z2(:,:) &
              - 2d0*(2d0*OVOV(:,:) - OVVO(:,:) - 2d0*OOVV(:,:)*t2(:,:))*z2(:,:)
 
      do i=1,O
        do a=1,V
 
          do j=1,O
            r2(i,a) = r2(i,a) - 2d0*OOVV(j,a)*t2(j,a)*z2(i,a) - 2d0*OOVV(i,a)*z2(j,a)*t2(j,a) & 
                              + OOOO(i,j)*z2(j,a) + yO(i,j)*z2(j,a)
          end do
 
          do b=1,V
            r2(i,a) = r2(i,a) - 2d0*OOVV(i,b)*t2(i,b)*z2(i,a) - 2d0*OOVV(i,a)*z2(i,b)*t2(i,b) & 
                              + VVVV(b,a)*z2(i,b) + yV(a,b)*z2(i,b)
          end do
 
        end do
      end do
 
     ! Check convergence 
 
      CvgAmp = maxval(abs(r2(:,:)))
    
     ! Update amplitudes
 
     z2(:,:) = z2(:,:) - 0.5d0*r2(:,:)/delta_OV(:,:)
 
     ! DIIS extrapolation
 
     if(max_diis > 1) then
 
       n_diis = min(n_diis+1,max_diis)
       call DIIS_extrapolation(rcond,O*V,O*V,n_diis,err_diis,z2_diis,-0.5d0*r2/delta_OV,z2)
 
      end if
 
      write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
        '|',nItAmp,'|',ECC+ENuc,'|',EcCC,'|',CvgAmp,'|'
 
    end do
    write(*,*)'----------------------------------------------------'
    write(*,*)

    !---------------------------!
    ! End Loop for z ampltiudes !
    !---------------------------!

    deallocate(r2,yO,yV)
    deallocate(err_diis,z2_diis)
 
    ! Did it actually converge?
 
    if(nItAmp == maxIt) then
 
      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'! Convergence failed for z ampltiudes !'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
 
      stop
  
    end if
 
    !--------------------------!
    ! Compute density matrices !
    !--------------------------!
 
    allocate(rdm1(N,N),rdm2(N,N,N,N))
    allocate(xOO(O,O),xVV(V,V),xOV(O,V))

    xOO(:,:) = matmul(t2,transpose(z2))
    xVV(:,:) = matmul(transpose(z2),t2)
    xOV(:,:) = matmul(t2,matmul(transpose(z2),t2))

    ! Form 1RDM
 
    rdm1(:,:) = 0d0
 
    do i=1,O
      rdm1(i,i) = 2d0*(1d0 - xOO(i,i))
    end do
 
    do a=1,V
      rdm1(O+a,O+a) = 2d0*xVV(a,a)
    end do
 
    ! Check 1RDM
 
    tr_1rdm = trace_matrix(N,rdm1)
    write(*,'(A25,F16.10)') ' --> Trace of the 1RDM = ',tr_1rdm
 
    if( abs(dble(2*O) - tr_1rdm) > thresh ) & 
    write(*,*) ' !!! Your 1RDM seems broken !!! '
    write(*,*)
 
    ! write(*,*) '1RDM is diagonal at the pCCD level:'
    ! call matout(N,N,rdm1)
 
    ! Form 2RM
 
    rdm2(:,:,:,:) = 0d0
 
    ! iijj
 
    do i=1,O
      do j=1,O
        rdm2(i,i,j,j) = 2d0*xOO(i,j)
      end do
    end do
 
    ! iiaa
 
    do i=1,O
      do a=1,V
        rdm2(i,i,O+a,O+a) = 2d0*(t2(i,a) + xOV(i,a) - 2d0*t2(i,a)*(xVV(a,a) + xOO(i,i) - t2(i,a)*z2(i,a)))
      end do
    end do
 
    ! aaii
 
    do i=1,O
      do a=1,V
        rdm2(O+a,O+a,i,i) = 2d0*z2(i,a)
      end do
    end do
 
    ! aabb
 
    do a=1,V
      do b=1,V
        rdm2(O+a,O+a,O+b,O+b) = 2d0*xVV(a,b)
      end do
    end do
 
    ! ijij
 
    do i=1,O
      do j=1,O
        rdm2(i,j,i,j) = 4d0*(1d0 - xOO(i,i) - xOO(j,j))
      end do
    end do
 
    ! ijji
 
    do i=1,O
      do j=1,O
        rdm2(i,j,j,i) = - 2d0*(1d0 - xOO(i,i) - xOO(j,j))
      end do
    end do
 
    ! iiii
 
    do i=1,O
      rdm2(i,i,i,i) = 2d0*(1d0 - xOO(i,i))
    end do
 
    ! iaia
 
    do i=1,O
      do a=1,V
        rdm2(i,O+a,i,O+a) = 4d0*(xVV(a,a) - t2(i,a)*z2(i,a))
      end do
    end do
 
    ! iaai
 
    do i=1,O
      do a=1,V
        rdm2(i,O+a,O+a,i) = - 2d0*(xVV(a,a) - t2(i,a)*z2(i,a))
      end do
    end do
 
    ! aiai
 
    do i=1,O
      do a=1,V
        rdm2(O+a,i,O+a,i) = 4d0*(xVV(a,a) - t2(i,a)*z2(i,a))
      end do
    end do
 
    ! aiia
 
    do i=1,O
      do a=1,V
        rdm2(O+a,i,i,O+a) = - 2d0*(xVV(a,a) - t2(i,a)*z2(i,a))
      end do
    end do
 
    ! abab
 
    do a=1,V
      rdm2(O+a,O+a,O+a,O+a) = 2d0*xVV(a,a)
    end do
 
    ! Check 2RDM
 
    tr_2rdm = trace_matrix(N**2,rdm2)
    write(*,'(A25,F16.10)') ' --> Trace of the 2RDM = ',tr_2rdm
 
    if( abs(dble(2*O*(2*O-1)) - tr_2rdm) > thresh ) & 
    write(*,*) ' !!! Your 2RDM seems broken !!! '
    write(*,*)
 
    ! write(*,*) '2RDM is not diagonal at the pCCD level:'
    ! call matout(N**2,N**2,rdm2)

    deallocate(xOO,xVV,xOV)
    deallocate(t2,z2)
 
    ! Compute electronic energy
 
    E1 = 0d0
    E2 = 0d0
 
    do p=1,N
      do q=1,N
        E1 = E1 + rdm1(p,q)*h(p,q)
        do r=1,N
          do s=1,N
            E2 = E2 + rdm2(p,q,r,s)*ERI_MO(p,q,r,s)
          end do
        end do
      end do
    end do
 
    E2 = 0.5d0*E2
 
    write(*,'(A25,F16.10)') ' One-electron energy = ',E1
    write(*,'(A25,F16.10)') ' Two-electron energy = ',E2
    write(*,'(A25,F16.10)') ' Electronic   energy = ',E1 + E2
    write(*,'(A25,F16.10)') ' Total        energy = ',E1 + E2 + ENuc
    write(*,*)
 
    !--------------------------!
    ! Compute orbital gradient !
    !--------------------------!
 
    allocate(grad(N**2))
 
    grad(:) = 0d0
 
    pq = 0
    do p=1,N
      do q=1,N
 
        pq = pq + 1
 
        do r=1,N
          grad(pq) = grad(pq) + h(r,p)*rdm1(r,q)  - h(q,r)*rdm1(p,r)
        end do
 
        do r=1,N
          do s=1,N
            do t=1,N
              grad(pq) = grad(pq) + (ERI_MO(r,s,p,t)*rdm2(r,s,q,t) - ERI_MO(q,t,r,s)*rdm2(p,t,r,s))
            end do
          end do
        end do
 
      end do
    end do
 
    write(*,*) 'Orbital gradient at the pCCD level:'
    call matout(N,N,grad)
    write(*,*)
 
   ! Check convergence of orbital optimization
 
    CvgOrb = maxval(abs(grad))
    write(*,*) ' Iteration',nItOrb,'for pCCD orbital optimization'
    write(*,*) ' Convergence of orbital gradient = ',CvgOrb
    write(*,*)
  
    !-------------------------!
    ! Compute orbital Hessian !
    !-------------------------!
 
    allocate(hess(N**2,N**2),tmp(N,N,N,N))
 
    tmp(:,:,:,:) = 0d0
 
    do p=1,N
      do q=1,N
 
        do r=1,N
          do s=1,N
 
            tmp(p,q,r,s) = - h(s,p)*rdm1(r,q) - h(q,r)*rdm1(p,s)
            
            do u=1,N
 
              tmp(p,q,r,s) = tmp(p,q,r,s) + 0.5d0*(                          &
                  Kronecker_delta(q,r)*(h(u,p)*rdm1(u,s) + h(s,u)*rdm1(p,u)) &
                + Kronecker_delta(p,s)*(h(u,r)*rdm1(u,q) + h(q,u)*rdm1(r,u)) )
 
            end do
 
            do u=1,N
              do w=1,N
 
              tmp(p,q,r,s) = tmp(p,q,r,s) + ERI_MO(u,w,p,r)*rdm2(u,w,q,s) + ERI_MO(q,s,u,w)*rdm2(p,r,u,w)
 
              end do
            end do
 
            do t=1,N
              do u=1,N
 
              tmp(p,q,r,s) = tmp(p,q,r,s) - (                                   &
                  ERI_MO(s,t,p,u)*rdm2(r,t,q,u) + ERI_MO(t,s,p,u)*rdm2(t,r,q,u) &
                + ERI_MO(q,u,r,t)*rdm2(p,u,s,t) + ERI_MO(q,u,t,r)*rdm2(p,u,t,s) )
 
              end do
            end do
 
            do t=1,N
              do u=1,N
                do w=1,N
 
                  tmp(p,q,r,s) = tmp(p,q,r,s) + 0.5d0*(                                                    & 
                      Kronecker_delta(q,r)*(ERI_MO(u,w,p,t)*rdm2(u,w,s,t) + ERI_MO(s,t,u,w)*rdm2(p,t,u,w)) &
                    + Kronecker_delta(p,s)*(ERI_MO(q,t,u,w)*rdm2(r,t,u,w) + ERI_MO(u,w,r,t)*rdm2(u,w,q,t)) )
 
                end do
              end do
            end do
 
          end do
        end do
 
      end do
    end do
 
    ! Flatten Hessian matrix and add permutations
 
    pq = 0
    do p=1,N
      do q=1,N
 
        pq = pq + 1
     
        rs = 0
        do r=1,N
          do s=1,N
 
            rs = rs + 1
 
            hess(pq,rs) = tmp(p,r,q,s) - tmp(r,p,q,s) - tmp(p,r,s,q) + tmp(r,p,s,q)
!!          hess(pq,rs) = tmp(p,q,r,s) - tmp(q,p,r,s) - tmp(p,q,s,r) + tmp(q,p,s,r)
 
          end do
        end do
 
      end do
    end do

    deallocate(rdm1,rdm2,tmp)
 
    allocate(hessInv(N**2,N**2))
 
    call inverse_matrix(N**2,hess,hessInv)

    deallocate(hess)
 
    allocate(Kap(N,N))
 
    Kap(:,:) = 0d0
 
    pq = 0
    do p=1,N
      do q=1,N
 
        pq = pq + 1
 
        rs = 0
        do r=1,N
          do s=1,N
 
            rs = rs + 1
 
              Kap(p,q) = Kap(p,q) - hessInv(pq,rs)*grad(rs)
 
          end do
        end do
 
      end do
    end do
 
    deallocate(hessInv,grad)

    write(*,*) 'kappa'
    call matout(N,N,Kap)
    write(*,*)
 
    allocate(ExpKap(N,N))
    call matrix_exponential(N,Kap,ExpKap)
    deallocate(Kap)
 
    write(*,*) 'e^kappa'
    call matout(N, N, ExpKap)
    write(*,*)
 
    write(*,*) 'Old orbitals'
    call matout(nBas, N, c)
    write(*,*)
 
    c = matmul(c, ExpKap)
    deallocate(ExpKap)
 
    write(*,*) 'Rotated orbitals'
    call matout(nBas, N, c)
    write(*,*)

  end do

  !-----------------------------------!
  ! End Loop for orbital optimization !
  !-----------------------------------!

    ! Did it actually converge?

    if(nItOrb == maxIt) then

      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'! Convergence failed for orbital optimization !'
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

      stop

    end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','pCCD correlation energy',EcCC)

  end if

end subroutine 
