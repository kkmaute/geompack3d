function minang ( xr, yr, xs, ys, ind, alpha, theta, vcl, pvl, iang )

!*****************************************************************************80
!
!! MINANG determines the minimum of the boundary angles for a separator.
!
!  Discussion:
!
!    This routine determines the minimum of the 4 angles at the boundary
!    resulting from using edge joining vertices (XR,YR) and
!    (XS,YS) as a separator.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XR, YR, the coordinates of the reflex vertex.
!
!    Input, real ( kind = 8 ) XS, YS, the coordinates of other endpoint of
!    possible separator.
!
!    Input, integer ( kind = 4 ) IND, if positive then (XS,YS) has index IND in PVL; else
!    (XS,YS) is on edge joining vertices with indices -IND
!    and SUCC(-IND) in PVL.
!
!    Input, real ( kind = 8 ) ALPHA, the polar angle of (XS,YS) with respect
!    to (XR,YR).
!
!    Input, real ( kind = 8 ) THETA, the interior angle at reflex vertex.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), real ( kind = 8 ) IANG(1:*), the polygon
!    vertex list, interior angles.
!
!    Output, real ( kind = 8 ) MINANG, the minimum of the 4 angles in radians.
!
  implicit none

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angle
  real    ( kind = 8 ) beta1
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  real    ( kind = 8 ) minang
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) theta
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) xr
  real    ( kind = 8 ) xs
  real    ( kind = 8 ) yr
  real    ( kind = 8 ) ys

  if ( 0 < ind ) then
    j = pvl(succ,ind)
    ang = iang(ind)
  else
    j = pvl(succ,-ind)
    ang = pi
  end if

  l = pvl(loc,j)
  beta1 = angle ( xr, yr, xs, ys, vcl(1,l), vcl(2,l) )

  minang = min ( alpha, theta - alpha, ang - beta1, beta1 )

  return
end
