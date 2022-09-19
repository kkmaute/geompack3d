function radrth ( a, b, c, d )

!*****************************************************************************80
!
!! RADRTH computes the aspect ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes the radius ratio or aspect ratio of a tetrahedron
!    = 3 * inradius / circumradius.
!
!  Modified:
!
!    16 August 2005
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), the vertices
!    of the tetrahedron.
!
!    Output, real ( kind = 8 ) RADRTH, the radius ratio of the tetrahedron.
!
  implicit none

  integer ( kind = 4 ), parameter :: dim_num = 3

  real    ( kind = 8 ) a(dim_num)
  real    ( kind = 8 ) ab(dim_num)
  real    ( kind = 8 ) ac(dim_num)
  real    ( kind = 8 ) ad(dim_num)
  real    ( kind = 8 ) b(dim_num)
  real    ( kind = 8 ) bc(dim_num)
  real    ( kind = 8 ) bd(dim_num)
  real    ( kind = 8 ) c(dim_num)
  real    ( kind = 8 ) cd(dim_num)
  real    ( kind = 8 ) cp1
  real    ( kind = 8 ) cp2
  real    ( kind = 8 ) cp3
  real    ( kind = 8 ) d(dim_num)
  real    ( kind = 8 ) denom
  real    ( kind = 8 ) fa
  real    ( kind = 8 ) fb
  real    ( kind = 8 ) fc
  real    ( kind = 8 ) fd
  real    ( kind = 8 ) lab
  real    ( kind = 8 ) lac
  real    ( kind = 8 ) lad
  real    ( kind = 8 ) lbc
  real    ( kind = 8 ) lbd
  real    ( kind = 8 ) lcd
  real    ( kind = 8 ) pb
  real    ( kind = 8 ) pc
  real    ( kind = 8 ) pd
  real    ( kind = 8 ) radrth
  real    ( kind = 8 ) t1
  real    ( kind = 8 ) t2
  real    ( kind = 8 ) vol
!
!  Compute the vectors that represent the sides.
!
  ab(1:dim_num) = b(1:dim_num) - a(1:dim_num)
  ac(1:dim_num) = c(1:dim_num) - a(1:dim_num)
  ad(1:dim_num) = d(1:dim_num) - a(1:dim_num)
  bc(1:dim_num) = c(1:dim_num) - b(1:dim_num)
  bd(1:dim_num) = d(1:dim_num) - b(1:dim_num)
  cd(1:dim_num) = d(1:dim_num) - c(1:dim_num)
!
!  Compute the squares of the lengths of the sides.
!
  lab = sum ( ab(1:dim_num)**2 )
  lac = sum ( ac(1:dim_num)**2 )
  lad = sum ( ad(1:dim_num)**2 )
  lbc = sum ( bc(1:dim_num)**2 )
  lbd = sum ( bd(1:dim_num)**2 )
  lcd = sum ( cd(1:dim_num)**2 )

  pb = sqrt ( lab * lcd )
  pc = sqrt ( lac * lbd )
  pd = sqrt ( lad * lbc )

  cp1 = ab(2) * ac(3) - ab(3) * ac(2)
  cp2 = ab(3) * ac(1) - ab(1) * ac(3)
  cp3 = ab(1) * ac(2) - ab(2) * ac(1)
  fd = sqrt ( cp1**2 + cp2**2 + cp3**2 )

  cp1 = ab(2) * ad(3) - ab(3) * ad(2)
  cp2 = ab(3) * ad(1) - ab(1) * ad(3)
  cp3 = ab(1) * ad(2) - ab(2) * ad(1)
  fc = sqrt ( cp1**2 + cp2**2 + cp3**2 )

  cp1 = bc(2) * bd(3) - bc(3) * bd(2)
  cp2 = bc(3) * bd(1) - bc(1) * bd(3)
  cp3 = bc(1) * bd(2) - bc(2) * bd(1)
  fa = sqrt ( cp1**2 + cp2**2 + cp3**2 )

  cp1 = ac(2) * ad(3) - ac(3) * ad(2)
  cp2 = ac(3) * ad(1) - ac(1) * ad(3)
  cp3 = ac(1) * ad(2) - ac(2) * ad(1)
  fb = sqrt ( cp1**2 + cp2**2 + cp3**2 )

  t1 = pb + pc
  t2 = pb - pc

  denom = ( fa + fb + fc + fd ) &
    * sqrt ( abs ( ( t1 + pd ) * ( t1 - pd ) * ( pd + t2 ) * ( pd - t2 ) ) )

  if ( denom == 0.0D+00 ) then
    radrth = 0.0D+00
  else
    vol = ab(1) * cp1 + ab(2) * cp2 + ab(3) * cp3
    radrth = 12.0D+00 * vol**2 / denom
  end if

  return
end
