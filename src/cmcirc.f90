function cmcirc ( x0, y0, x1, y1, x2, y2, x3, y3 )

!*****************************************************************************80
!
!! CMCIRC determines if a point is in the circumcircle of three points.
!
!  Discussion:
!
!    This routine determines whether (X0,Y0) is in the circumcircle through
!    the three points (X1,Y1), (X2,Y2), (X3,Y3).
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
!    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
!    vertex coordinates.
!
!    Output, integer ( kind = 4 ) CMCIRC:
!     2 if three vertices collinear,
!     1 if ( X0,Y0) inside circle,
!     0 if ( X0,Y0) on circle,
!    -1 if ( X0,Y0) outside circle
!
  implicit none

  real    ( kind = 8 ) a11
  real    ( kind = 8 ) a12
  real    ( kind = 8 ) a21
  real    ( kind = 8 ) a22
  real    ( kind = 8 ) b1
  real    ( kind = 8 ) b2
  integer ( kind = 4 ) cmcirc
  real    ( kind = 8 ) det
  real    ( kind = 8 ) diff
  real    ( kind = 8 ) rsq
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) x2
  real    ( kind = 8 ) x3
  real    ( kind = 8 ) xc
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) y2
  real    ( kind = 8 ) y3
  real    ( kind = 8 ) yc
!
  tol = 100.0D+00 * epsilon ( tol )
  cmcirc = 2
  a11 = x2 - x1
  a12 = y2 - y1
  a21 = x3 - x1
  a22 = y3 - y1
  tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = a11**2 + a12**2
  b2 = a21**2 + a22**2
  det = det + det
  xc = ( b1 * a22 - b2 * a12 ) / det
  yc = ( b2 * a11 - b1 * a21 ) / det
  rsq = xc**2 + yc**2
  diff = (( x0 - x1 - xc )**2 + ( y0 - y1 - yc)**2 ) - rsq
  tolabs = tol*rsq

  if ( diff < -tolabs ) then
    cmcirc = 1
  else if ( tolabs < diff ) then
    cmcirc = -1
  else
    cmcirc = 0
  end if

  return
end
