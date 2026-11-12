subroutine TDA_Faddeev_static3(nERI,nOrb,nphys,nC,nO,nR,eHF,ERI,sigma)
! Strict MP third-order static self-energy in the canonical HF basis.
! With spatial ERIs, use interleaved spin orbitals and return the alpha block.
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer,intent(in) :: nERI,nOrb,nphys,nC,nO,nR
  double precision,intent(in) :: eHF(nOrb),ERI(nERI,nERI,nERI,nERI)
  double precision,intent(out) :: sigma(nphys,nphys)
  double precision,allocatable :: t(:,:,:,:),density(:,:)
  double precision :: value,den
  integer :: i,j,k,a,b,c,p,q,r,s,step,first,last

  if(nOrb/=nERI.and.nOrb/=2*nERI) error stop 'Static self-energy: invalid ERI dimension'
  if(nOrb/=nphys.and.nOrb/=2*nphys) error stop 'Static self-energy: invalid physical dimension'
  first=nC+1; last=nOrb-nR; step=nOrb/nphys
  allocate(t(first:nO,first:nO,nO+1:last,nO+1:last),density(first:last,first:last))

  do b=nO+1,last
    do a=nO+1,last
      do j=first,nO
        do i=first,nO
          den=eHF(i)+eHF(j)-eHF(a)-eHF(b)
          if(abs(den)<1d-12) error stop 'Static self-energy: singular doubles denominator'
          t(i,j,a,b)=v(a,b,i,j)/den
        end do
      end do
    end do
  end do

  do j=first,nO
    do i=first,j
      density(i,j)=-0.5d0*sum(t(i,:,:,:)*t(j,:,:,:))
      density(j,i)=density(i,j)
    end do
  end do
  do b=nO+1,last
    do a=nO+1,b
      density(a,b)=0.5d0*sum(t(:,:,a,:)*t(:,:,b,:))
      density(b,a)=density(a,b)
    end do
  end do

  do a=nO+1,last
    do i=first,nO
      value=0d0
      do c=nO+1,last
        do b=nO+1,last
          do j=first,nO
            value=value+v(a,j,b,c)*t(i,j,b,c)
          end do
        end do
      end do
      do b=nO+1,last
        do k=first,nO
          do j=first,nO
            value=value-v(j,k,i,b)*t(j,k,a,b)
          end do
        end do
      end do
      den=eHF(i)-eHF(a)
      if(abs(den)<1d-12) error stop 'Static self-energy: singular singles denominator'
      density(a,i)=0.5d0*value/den
      density(i,a)=density(a,i)
    end do
  end do

  do q=1,nphys
    do p=1,q
      value=0d0
      do s=first,last
        do r=first,last
          value=value+v(1+step*(p-1),r,1+step*(q-1),s)*density(s,r)
        end do
      end do
      sigma(p,q)=value
      sigma(q,p)=value
    end do
  end do
  if(.not.all(ieee_is_finite(sigma))) error stop 'Nonfinite third-order static self-energy'
  write(*,*) 'FTDA: third-order static self-energy enabled'

contains

  double precision function v(p,q,r,s)
    integer,intent(in) :: p,q,r,s
    if(nERI==nOrb) then
      v=ERI(p,q,r,s)-ERI(p,q,s,r)
    else
      v=0d0
      if(mod(p-r,2)==0.and.mod(q-s,2)==0) v=ERI((p+1)/2,(q+1)/2,(r+1)/2,(s+1)/2)
      if(mod(p-s,2)==0.and.mod(q-r,2)==0) v=v-ERI((p+1)/2,(q+1)/2,(s+1)/2,(r+1)/2)
    end if
  end function
end subroutine
