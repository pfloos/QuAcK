subroutine qsGW_PT(nBas,nC,nO,nV,nR,nS,e0,SigCm)

! Compute the 1st-, 2nd-, 3rd- and 4th-order correction on the qsGW quasiparticle energies

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e0(nBas),SigCm(nBas,nBas)

! Local variables

  integer                       :: x,y,z,t
  double precision              :: eps
  double precision,allocatable  :: e1(:),e2(:),e3(:),e4(:)
  double precision,parameter    :: threshold = 1d-15

! Allocation

  allocate(e1(nBas),e2(nBas),e3(nBas),e4(nBas))

! Initalization

  e1(:) = 0d0
  e2(:) = 0d0
  e3(:) = 0d0
  e4(:) = 0d0

! Print zeroth-order qsGW QP energies 

  write(*,*)
  write(*,'(A50)') '-----------------------------------------------'
  write(*,'(A50)') ' 0th-order values of qsGW QP energies (eV) '
  write(*,'(A50)') '-----------------------------------------------'
  call matout(nBas,1,e0(:)*HaToeV)

! Compute 1st-order correction of qsGW QP energies 

  do x=nC+1,nBas-nR

    e1(x) = SigCm(x,x)

  end do

  write(*,*)
  write(*,'(A50)') '-----------------------------------------------'
  write(*,'(A50)') ' 1st-order correction of qsGW QP energies (eV) '
  write(*,'(A50)') '-----------------------------------------------'
  call matout(nBas,1,e1(:)*HaToeV)

! Compute 2nd-order correction of qsGW QP energies 

  do x=nC+1,nBas-nR
    do y=nC+1,nBas-nR

      eps = e0(x) - e0(y)
      if(abs(eps) > threshold) e2(x) = e2(x) + SigCm(x,y)**2/eps

    end do
  end do

  write(*,*)
  write(*,'(A50)') '-----------------------------------------------'
  write(*,'(A50)') ' 2nd-order correction of qsGW QP energies (eV) '
  write(*,'(A50)') '-----------------------------------------------'
  call matout(nBas,1,e2(:)*HaToeV)

! Compute 3nd-order correction of qsGW QP energies 

  do x=nC+1,nBas-nR
    do y=nC+1,nBas-nR
      do z=nC+1,nBas-nR

        eps = (e0(x) - e0(y))*(e0(x) - e0(z))
        if(abs(eps) > threshold) e3(x) = e3(x) + SigCm(x,y)*SigCm(y,z)*SigCm(z,x)/eps

      end do
    end do
  end do

  write(*,*)
  write(*,'(A50)') '-----------------------------------------------'
  write(*,'(A50)') ' 3rd-order correction of qsGW QP energies (eV) '
  write(*,'(A50)') '-----------------------------------------------'
  call matout(nBas,1,e3(:)*HaToeV)

! Compute 4nd-order correction of qsGW QP energies 

  do x=nC+1,nBas-nR
    do y=nC+1,nBas-nR
      do z=nC+1,nBas-nR
        do t=nC+1,nBas-nR

          eps = (e0(x) - e0(y))*(e0(x) - e0(z))*(e0(x) - e0(t))
          if(abs(eps) > threshold) e4(x) = e4(x) + SigCm(x,y)*SigCm(y,z)*SigCm(z,t)*SigCm(t,x)/eps

        end do
      end do
    end do
  end do

  do x=nC+1,nBas-nR
    do y=nC+1,nBas-nR

      eps = (e0(x) - e0(y))**2
      if(abs(eps) > threshold) e4(x) = e4(x) - e2(x)*SigCm(x,y)**2/eps

    end do
  end do

  write(*,*)
  write(*,'(A50)') '-----------------------------------------------'
  write(*,'(A50)') ' 4th-order correction of qsGW QP energies (eV) '
  write(*,'(A50)') '-----------------------------------------------'
  call matout(nBas,1,e4(:)*HaToeV)

end subroutine qsGW_PT
