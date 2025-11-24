subroutine Sigmac_MO_RHF_GF2_analytical(nBas,nOrb,nO,verbose,c,eHF,nfreqs,wcoord,ERI_AO,Sigma_c,Ec)


! Use the restricted Sigma_c(E) to compute the linnearized approximation to G

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: verbose
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO

  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(in)   :: wcoord(nfreqs)
  double precision,intent(in)   :: c(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ifreq
  integer                       :: iorb,jorb,aorb,borb,porb,qorb

  complex*16                    :: num,eps

  double precision,allocatable  :: ERI_MO(:,:,:,:)

! Ouput variables

  double precision,intent(out)  :: Ec
  complex *16,intent(out)       :: Sigma_c(nfreqs,nOrb,nOrb)

! Allocate and initialize arrays and variables
  allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
  call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

  Sigma_c=czero

! Sigma_c

  do ifreq=1,nfreqs

    do porb=1,nOrb
      do qorb=1,nOrb
        do iorb=1,nO
          do jorb=1,nO
            do aorb=nO+1,nOrb
              eps = im*wcoord(ifreq) + eHF(aorb) - eHF(iorb) - eHF(jorb)
              num = (2d0*ERI_MO(porb,aorb,iorb,jorb) - ERI_MO(porb,aorb,jorb,iorb))*ERI_MO(qorb,aorb,iorb,jorb)
              Sigma_c(ifreq,porb,qorb) = Sigma_c(ifreq,porb,qorb) + num/eps
            end do
          end do
        end do
      end do
    end do

    do porb=1,nOrb
      do qorb=1,nOrb
        do iorb=1,nO
          do aorb=nO+1,nOrb
            do borb=nO+1,nOrb
              eps = im*wcoord(ifreq) + eHF(iorb) - eHF(aorb) - eHF(borb)
              num = (2d0*ERI_MO(porb,iorb,aorb,borb) - ERI_MO(porb,iorb,borb,aorb))*ERI_MO(qorb,iorb,aorb,borb)
              Sigma_c(ifreq,porb,qorb) = Sigma_c(ifreq,porb,qorb) + num/eps
            end do
          end do
        end do
      end do
    end do

   if(verbose/=0) then
    write(*,'(a,*(f20.8))') ' Analytic ',im*wcoord(ifreq)
    do iorb=1,nOrb
     write(*,'(*(f20.8))') Sigma_c(ifreq,iorb,:)
    enddo
   endif
  
  enddo

  Ec = 0d0
  do iorb=1,nO
    do jorb=1,nO
      do aorb=nO+1,nOrb
        do borb=nO+1,nOrb
          eps = eHF(iorb) + eHF(jorb) - eHF(aorb) - eHF(borb)
          num = (2d0*ERI_MO(iorb,jorb,aorb,borb) - ERI_MO(iorb,jorb,borb,aorb))*ERI_MO(iorb,jorb,aorb,borb)
          Ec = Ec + real(num/eps)
        end do
      end do
    end do
  end do

  deallocate(ERI_MO)

end subroutine

