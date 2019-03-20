subroutine S3eInt(debug,iType,np3eInt,nSigp3eInt, &
                  ExpS,KG,DG,ExpG,                    &
                  ExpBra,CenterBra,AngMomBra,         &
                  ExpKet,CenterKet,AngMomKet,         &
                  p3eInt)

! Perform contraction over the operators for three-electron integrals

  implicit none
  include 'parameters.h'


! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: iType
  double precision,intent(in)   :: ExpS
  integer,intent(in)            :: KG
  double precision,intent(in)   :: DG(KG),ExpG(KG)
  double precision,intent(in)   :: ExpBra(3),ExpKet(3)
  double precision,intent(in)   :: CenterBra(3,3),CenterKet(3,3)
  integer,intent(in)            :: AngMomBra(3,3),AngMomKet(3,3)

! Local variables

  double precision              :: ExpSG13,ExpSG23
  double precision              :: G3eInt

  integer                       :: k,l

! Output variables

  integer,intent(out)           :: np3eInt,nSigp3eInt
  double precision              :: p3eInt

  p3eInt = 0d0

  if(iType == 1) then
  
    do k=1,KG
      ExpSG13 = ExpG(k)*ExpS**2
      p3eInt = p3eInt                                   &
             + DG(k)*G3eInt(debug,iType,                &
                            ExpSG13,ExpSG23,            &
                            ExpBra,CenterBra,AngMomBra, &
                            ExpKet,CenterKet,AngMomKet)
    end do

  end if

  if(iType == 2 .or. iType == 3) then
  
    do k=1,KG
      do l=1,KG
        ExpSG13 = ExpG(k)*ExpS**2
        ExpSG23 = ExpG(l)*ExpS**2
        p3eInt = p3eInt                                         &
               + DG(k)*DG(l)*G3eInt(debug,iType,                &
                                    ExpSG13,ExpSG23,            &
                                    ExpBra,CenterBra,AngMomBra, &
                                    ExpKet,CenterKet,AngMomKet)
      end do
    end do

  end if

! Print result

  np3eInt = np3eInt + 1

  if(abs(p3eInt) > 1d-15) then
    nSigp3eInt = nSigp3eInt + 1
    if(.false.) write(*,'(A15,1X,F16.10)') '[a1a2a3|b1b2b3] = ',p3eInt
  end if

end subroutine S3eInt
