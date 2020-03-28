  integer,parameter              :: ncart = 3
  integer,parameter              :: nspin = 2
  integer,parameter              :: nsp   = 3
  integer,parameter              :: maxEns = 10
  integer,parameter              :: maxShell = 512
  integer,parameter              :: maxL = 7
  integer,parameter              :: n1eInt = 3
  integer,parameter              :: n2eInt = 4
  integer,parameter              :: n3eInt = 3
  integer,parameter              :: n4eInt = 3
  integer,parameter              :: maxK = 20

  double precision,parameter    :: threshold = 1d-15
  double precision,parameter    :: pi = acos(-1d0)
  double precision,parameter    :: HaToeV = 27.21138602d0
  double precision,parameter    :: pmtoau = 0.0188973d0
  double precision,parameter    :: BoToAn = 0.529177249d0

  double precision,parameter    :: CxLDA = - (3d0/4d0)*(3d0/pi)**(1d0/3d0)
  double precision,parameter    :: Cx0   = - (4d0/3d0)*(1d0/pi)**(1d0/3d0)
! double precision,parameter    :: Cx1   = - 0.904d0*(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  double precision,parameter    :: Cx1   = - (176d0/105d0)*(1d0/pi)**(1d0/3d0)

