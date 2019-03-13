function element_number(element_name)

 implicit none

 integer,parameter                  :: nelement_max = 103
 character(len=2),intent(in)        :: element_name
 integer                            :: element_number
 character(len=2),parameter         :: element_list(nelement_max) =                            &
  (/' H',                                                                                'He', &  !   2
    'Li','Be',                                                  ' B',' C',' N',' O',' F','Ne', &  !  10
    'Na','Mg',                                                  'Al','Si',' P',' S','Cl','Ar', &  !  18
    ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  !  36
    'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe', &  !  54
    'Cs','Ba',                                                                                 &  !  56
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',                     &  !  70
              'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &  !  86
    'Fr','Ra',                                                                                 &  !  88
              'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',           &  ! 102
    'Lr'                                                                                       &  ! 103
  /)

!=====
 integer :: ielement
!=====

 ielement=1
 do while( ADJUSTL(element_name) /= ADJUSTL(element_list(ielement)) )
   if( ielement == nelement_max ) then
     write(*,'(a,a)')    ' Input symbol ',element_name
     write(*,'(a,i3,a)') ' Element symbol is not one of first ',nelement_max,' elements'
     write(*,*) '!!! element symbol not understood !!!'
     stop
   endif
   ielement = ielement + 1
 enddo

 element_number = ielement

end function element_number

function element_core(zval,zatom)
 implicit none
 double precision,intent(in) :: zval
 double precision,intent(in) :: zatom
 integer                     :: element_core
!=====

 !
 ! If zval /= zatom, this is certainly an effective core potential
 ! and no core states should be frozen.
 if( ABS(zval - zatom) > 1d0-3 ) then
   element_core = 0
 else

   if( zval <= 4.00001d0 ) then  ! up to Be
     element_core = 0
   else if( zval <= 12.00001d0 ) then  ! up to Mg
     element_core = 1
   else if( zval <= 30.00001d0 ) then  ! up to Ca
     element_core = 5
   else if( zval <= 48.00001d0 ) then  ! up to Sr
     element_core = 9
   else
     write(*,*) '!!! not imlemented in element_core !!!'
     stop
   endif

 endif


end function element_core

function element_covalent_radius(zatom)

! Return covalent radius of an atom

 implicit none
 include 'parameters.h'

 integer,intent(in) :: zatom
 double precision   :: element_covalent_radius

 !
 ! Data from Cambridge Structural Database
 ! http://en.wikipedia.org/wiki/Covalent_radius
 !
 ! Values are first given in picometer
 ! They will be converted in bohr later on
 select case(zatom)
 case( 1)
   element_covalent_radius =  31.
 case( 2)
   element_covalent_radius =  28.
 case( 3)
   element_covalent_radius = 128.
 case( 4)
   element_covalent_radius =  96.
 case( 5)
   element_covalent_radius =  84.
 case( 6)
   element_covalent_radius =  73.
 case( 7)
   element_covalent_radius =  71.
 case( 8)
   element_covalent_radius =  66.
 case( 9)
   element_covalent_radius =  57.
 case(10) ! Ne.
   element_covalent_radius =  58.
 case(11)
   element_covalent_radius = 166.
 case(12)
   element_covalent_radius = 141.
 case(13)
   element_covalent_radius = 121.
 case(14)
   element_covalent_radius = 111.
 case(15)
   element_covalent_radius = 107.
 case(16)
   element_covalent_radius = 105.
 case(17)
   element_covalent_radius = 102.
 case(18) ! Ar.
   element_covalent_radius = 106.
 case(19)
   element_covalent_radius = 203.
 case(20)
   element_covalent_radius = 176.
 case(21)
   element_covalent_radius = 170.
 case(22)
   element_covalent_radius = 160.
 case(23)
   element_covalent_radius = 153.
 case(24)
   element_covalent_radius = 139.
 case(25)
   element_covalent_radius = 145.
 case(26)
   element_covalent_radius = 145.
 case(27)
   element_covalent_radius = 140.
 case(28)
   element_covalent_radius = 124.
 case(29)
   element_covalent_radius = 132.
 case(30)
   element_covalent_radius = 122.
 case(31)
   element_covalent_radius = 120.
 case(32)
   element_covalent_radius = 119.
 case(34)
   element_covalent_radius = 120.
 case(35)
   element_covalent_radius = 120.
 case(36) ! Kr.
   element_covalent_radius = 116.
 case default
   write(*,*) '!!! covalent radius not available !!!'
   stop
 end select

 ! pm to bohr conversion
 element_covalent_radius = element_covalent_radius*pmtoau


end function element_covalent_radius

