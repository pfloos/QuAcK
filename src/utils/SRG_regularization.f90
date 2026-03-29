double precision function SRG_reg(dem,flow)

  implicit none

  double precision,intent(in)  :: dem,flow

  SRG_reg = (1d0 - exp(-2d0*flow*dem*dem))/dem

end function 

double precision function SRG_reg2(dem1,dem2,flow)

  implicit none

  double precision,intent(in)  :: dem1,dem2,flow

  SRG_reg2 = (1d0 - exp(-2d0*flow*(dem1+dem2)**2))/(dem1*dem2)

end function 
