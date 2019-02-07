subroutine S2eInt(debug,iType,np2eInt,nSigp2eInt, &
                  ExpS,KG,DG,ExpG,                &
                  ExpBra,CenterBra,AngMomBra,     &
                  ExpKet,CenterKet,AngMomKet,     &
                  p2eInt)

! Perform contraction over the operator for two-electron integrals

  implicit none
  include 'parameters.h'


! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: iType
  double precision,intent(in)   :: ExpS
  integer,intent(in)            :: KG
  double precision,intent(in)   :: DG(KG),ExpG(KG)
  double precision,intent(in)   :: ExpBra(2),ExpKet(2)
  double precision,intent(in)   :: CenterBra(2,3),CenterKet(2,3)
  integer,intent(in)            :: AngMomBra(2,3),AngMomKet(2,3)

! Local variables

  double precision              :: ExpSG
  double precision              :: G2eInt,GF12Int

  integer                       :: k

! Output variables

  integer,intent(out)           :: np2eInt,nSigp2eInt
  double precision              :: p2eInt

  p2eInt = 0d0

! Gaussian geminal

  if(iType == 2) then
    do k=1,KG
      ExpSG = ExpG(k)*ExpS**2
      p2eInt = p2eInt                                                     &
             + DG(k)*GF12Int(ExpSG,                                       &
                             ExpBra(1),CenterBra(1,1:3),AngMomBra(1,1:3), &
                             ExpKet(1),CenterKet(1,1:3),AngMomKet(1,1:3), &
                             ExpBra(2),CenterBra(2,1:3),AngMomBra(2,1:3), &
                             ExpKet(2),CenterKet(2,1:3),AngMomKet(2,1:3))
    enddo
  else
    do k=1,KG
      ExpSG = ExpG(k)*ExpS**2
      p2eInt = p2eInt                                   &
             + DG(k)*G2eInt(debug,iType,                &
                            ExpSG,                      &
                            ExpBra,CenterBra,AngMomBra, &
                            ExpKet,CenterKet,AngMomKet)
    enddo
  endif

! Print result

  np2eInt = np2eInt + 1

  if(abs(p2eInt) > 1d-15) then
    nSigp2eInt = nSigp2eInt + 1
    if(.false.) write(*,'(A15,1X,F16.10)') '[a1a2|b1b2] = ',p2eInt
  endif

end subroutine S2eInt
