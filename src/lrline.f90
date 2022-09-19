function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

!*****************************************************************************80
!
!! LRLINE determines whether a point is left, right, or on a directed line.
!
!  Discussion:
!
!    This routine determines whether a point is to the left of, right of,
!    or on a directed line parallel to a line through given points.
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
!    Input, real ( kind = 8 ) XU, YU, XV1, YV1, XV2, YV2, vertex coordinates;
!    the directed line is parallel to and at signed distance DV to the
!    left of the directed line from (XV1,YV1) to (XV2,YV2);
!    (XU,YU) is the vertex for which the position
!    relative to the directed line is to be determined.
!
!    Input, real ( kind = 8 ) DV, signed distance (positive for left).
!
!    Output, integer ( kind = 4 ) LRLINE, +1, 0, or -1 depending on whether (XU,YU) is
!    to the right of, on, or left of the directed line
!    (0 if line degenerates to a point).
!
  implicit none

  real    ( kind = 8 ) dv
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dxu
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dyu
  integer ( kind = 4 ) lrline
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xv1
  real    ( kind = 8 ) xv2
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yv1
  real    ( kind = 8 ) yv2

  tol = 100.0D+00 * epsilon ( tol )
  dx = xv2 - xv1
  dy = yv2 - yv1
  dxu = xu - xv1
  dyu = yu - yv1

  tolabs = tol * max ( &
    abs ( dx ), abs ( dy ), abs ( dxu ), abs ( dyu ), abs ( dv ) )

  t = dy * dxu - dx * dyu

  if ( dv /= 0.0D+00 ) then
    t = t + dv * sqrt ( dx**2 + dy**2 )
  end if

  lrline = int ( sign ( 1.0D+00, t ) )

  if ( abs ( t ) <= tolabs ) then
    lrline = 0
  end if

  return
end
