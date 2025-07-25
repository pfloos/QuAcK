subroutine R_Parquet_QP_graph(eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,ERI, &
                              eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,   &
                              ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,   &
                              hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,   &
                              eHF,eQPlin,eOld,eQP,Z)

! Compute the graphical solution of the QP equation

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb,nC,nO,nV,nR
  integer,intent(in)            :: nS,nOOs,nVVs,nOOt,nVVt
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: eh_sing_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_sing_Om(nS)
  double precision,intent(in)   :: eh_trip_rho(nOrb,nOrb,nS)
  double precision,intent(in)   :: eh_trip_Om(nS)
  double precision,intent(in)   :: ee_sing_rho(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: ee_sing_Om(nVVs)
  double precision,intent(in)   :: ee_trip_rho(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: ee_trip_Om(nVVt)
  double precision,intent(in)   :: hh_sing_rho(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: hh_sing_Om(nOOs)
  double precision,intent(in)   :: hh_trip_rho(nOrb,nOrb,nOOt)
  double precision,intent(in)   :: hh_trip_Om(nOOt)

  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: eQPlin(nOrb)
  double precision,intent(in)   :: eOld(nOrb)

! Local variables

  integer                       :: p
  integer                       :: nIt
  integer,parameter             :: maxIt = 64
  double precision,parameter    :: thresh = 1d-6
  double precision              :: SigC,dSigC
  double precision              :: f,df
  double precision              :: w

! Output variables

  double precision,intent(out)  :: eQP(nOrb)
  double precision,intent(out)  :: Z(nOrb)

! Run Newton's algorithm to find the root

  write(*,*)'-----------------------------------------------------'
  write(*,'(A5,1X,A3,1X,A15,1X,A15,1X,A10)') 'Orb.','It.','e_QPlin (eV)','e_QP (eV)','Z'
  write(*,*)'-----------------------------------------------------'

  do p=nC+1,nOrb-nR

    w = eQPlin(p)
    nIt = 0
    f = 1d0
    
    do while (abs(f) > thresh .and. nIt < maxIt)
    
      nIt = nIt + 1

      call R_Parquet_self_energy_omega(p,w,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,ERI, &
                                       eh_sing_rho,eh_sing_Om,eh_trip_rho,eh_trip_Om,       &
                                       ee_sing_rho,ee_sing_Om,ee_trip_rho,ee_trip_Om,       &
                                       hh_sing_rho,hh_sing_Om,hh_trip_rho,hh_trip_Om,       &
                                       eOld,SigC,dSigC)

      f  = w - eHF(p) - SigC
      df = 1d0/(1d0 - dSigC)
      w = w - df*f
    
    end do
 
    if(nIt == maxIt) then 

      eQP(p) = eQPlin(p)
      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6,1X,A12)') p,nIt,eQPlin(p)*HaToeV,eQP(p)*HaToeV,Z(p),'Cvg Failed!'

    else

      eQP(p) = w
      Z(p)   = df

      write(*,'(I5,1X,I3,1X,F15.9,1X,F15.9,1X,F10.6)') p,nIt,eQPlin(p)*HaToeV,eQP(p)*HaToeV,Z(p)

    end if
          
  end do

  write(*,*)'-----------------------------------------------------'
  write(*,*)

end subroutine 
