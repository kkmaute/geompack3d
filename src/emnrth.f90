function emnrth ( a, b, c, d )

!*****************************************************************************80
!
!! EMNRTH computes the mean ratio of a tetrahedron.
!
!  Discussion:
!
!    This routine computes the eigenvalue or mean ratio of a tetrahedron
!    = 12 * ( 3 * volume )**(2/3) / (sum of square of edge lengths).
!
!    This value may be used as a shape quality measure for the tetrahedron.
!
!    For an equilateral tetrahedron, the value of this quality measure
!    will be 1.  For any other tetrahedron, the value will be between
!    0 and 1.
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
!    Output, real ( kind = 8 ) EMNRTH, the mean ratio of the tetrahedron.
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
  real    ( kind = 8 ) d(dim_num)
  real    ( kind = 8 ) denom
  real    ( kind = 8 ) emnrth
  real    ( kind = 8 ) lab
  real    ( kind = 8 ) lac
  real    ( kind = 8 ) lad
  real    ( kind = 8 ) lbc
  real    ( kind = 8 ) lbd
  real    ( kind = 8 ) lcd
  real    ( kind = 8 ) vol
!
!  Compute the vectors representing the sides of the tetrahedron.
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
!
!  Compute the volume.
!
  vol = abs ( &
      ab(1) * ( ac(2) * ad(3) - ac(3) * ad(2) ) &
    + ab(2) * ( ac(3) * ad(1) - ac(1) * ad(3) ) &
    + ab(3) * ( ac(1) * ad(2) - ac(2) * ad(1) ) ) / 6.0D+00

  denom = lab + lac + lad + lbc + lbd + lcd

  if ( denom == 0.0D+00 ) then
    emnrth = 0.0D+00
  else
    emnrth = 12.0D+00 * ( 3.0D+00 * vol )**( 2.0D+00 / 3.0D+00 ) / denom
  end if

  return
end
