double precision function SRG_reg(dem,flow)

  implicit none

  double precision,intent(in)  :: dem,flow

  SRG_reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

end function 

double precision function SRG_reg2(dem1,dem2,flow)

  implicit none

  double precision,intent(in)  :: dem1,dem2,flow

! SRG_reg2 = (1d0 - exp(-2d0*flow*(dem1+dem2)**2))/(dem1*dem2)
  SRG_reg2 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
  SRG_reg2 = SRG_reg2*(1d0 - exp(-2d0*flow*dem2*dem2))/dem2

end function 

double precision function SRG_reg3(dem1,dem2,dem3,flow)

  implicit none

  double precision,intent(in)  :: dem1,dem2,dem3,flow

  SRG_reg3 = (1d0 - exp(-2d0*flow*dem1*dem1))/dem1
  SRG_reg3 = SRG_reg3*(1d0 - exp(-2d0*flow*dem2*dem2))/dem2
  SRG_reg3 = SRG_reg3*(1d0 - exp(-2d0*flow*dem3*dem3))/dem3

end function 
