function angle ( xa, ya, xb, yb, xc, yc )

!*****************************************************************************80
!
!! ANGLE computes the size of an angle in 2D.
!
!  Discussion:
!
!    This routine computes the interior angle in radians at vertex
!    (XB,YB) of the chain formed by the directed edges from
!    (XA,YA) to (XB,YB) to (XC,YC).  The interior is to the
!    left of the two directed edges.
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
!    using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XA, YA, XB, YB, XC, YC, the vertex coordinates.
!
!    Output, real ( kind = 8 ) ANGLE, the angle, between 0 and 2*PI.
!    ANGLE is set to PI/2 in the degenerate case.
!
  implicit none

  real    ( kind = 8 ) angle
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) xa
  real    ( kind = 8 ) xb
  real    ( kind = 8 ) xc
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) ya
  real    ( kind = 8 ) yb
  real    ( kind = 8 ) yc

  tol = 100.0D+00 * epsilon ( tol )
  x1 = xa - xb
  y1 = ya - yb
  x2 = xc - xb
  y2 = yc - yb

  t = sqrt ( ( x1 * x1 + y1 * y1 ) * ( x2 * x2 + y2 * y2 ) )
  if ( t == 0.0D+00 ) then
    t = 1.0D+00
  end if

  t = ( x1 * x2 + y1 * y2 ) / t

  if ( 1.0D+00 - tol < abs ( t ) ) then
    t = sign ( 1.0D+00, t )
  end if

  angle = acos ( t )

  if ( x2 * y1 - y2 * x1 < 0.0D+00 ) then
    angle = 2.0D+00 * pi - angle
  end if

  return
end
