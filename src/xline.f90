subroutine xline ( xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, dv, dw, &
  xu, yu, parall )

!*****************************************************************************80
!
!! XLINE intersects two lines parallel to given lines.
!
!  Discussion:
!
!    This routine determines the intersection point of two lines parallel
!    to lines through given points.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2,
!    the vertex coordinates; first line is parallel to and at signed
!    distance DV to left of directed line from (XV1,YV1) to (XV2,YV2);
!    second line is parallel to and at signed distance DW to
!    left of directed line from (XW1,YW1) to (XW2,YW2).
!
!    Input, real ( kind = 8 ) DV, DW, the signed distances (positive for
!    left).
!
!    Output, real ( kind = 8 ) XU, YU, the coordinates of the point of
!    intersection iff PARALL is FALSE.
!
!    Output, logical PARALL, TRUE if the lines are parallel or two points for a
!    line are identical, FALSE otherwise.
!
  implicit none

  real    ( kind = 8 ) a11
  real    ( kind = 8 ) a12
  real    ( kind = 8 ) a21
  real    ( kind = 8 ) a22
  real    ( kind = 8 ) b1
  real    ( kind = 8 ) b2
  real    ( kind = 8 ) det
  real    ( kind = 8 ) dv
  real    ( kind = 8 ) dw
  logical              parall
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xv1
  real    ( kind = 8 ) xv2
  real    ( kind = 8 ) xw1
  real    ( kind = 8 ) xw2
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yv1
  real    ( kind = 8 ) yv2
  real    ( kind = 8 ) yw1
  real    ( kind = 8 ) yw2

  tol = 100.0D+00 * epsilon ( tol )
  parall = .true.
  a11 = yv2 - yv1
  a12 = xv1 - xv2
  a21 = yw2 - yw1
  a22 = xw1 - xw2
  tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
  det = a11 * a22 - a21 * a12

  if ( abs ( det ) <= tolabs ) then
    return
  end if

  b1 = xv1 * a11 + yv1 * a12

  if ( dv /= 0.0D+00 ) then
    b1 = b1 - dv * sqrt ( a11**2 + a12**2 )
  end if

  b2 = xw1 * a21 + yw1 * a22

  if ( dw /= 0.0D+00 ) then
    b2 = b2 - dw * sqrt ( a21**2 + a22**2 )
  end if

  xu = ( b1 * a22 - b2 * a12 ) / det
  yu = ( b2 * a11 - b1 * a21 ) / det
  parall = .false.

  return
end
