subroutine xedge ( mode, xv1, yv1, xv2, yv2, xw1, yw1, xw2, yw2, xu, yu, &
  intsct )

!*****************************************************************************80
!
!! XEDGE determines whether two edges, or an edge and a ray intersect.
!
!  Discussion:
!
!    This routine determines whether two edges or a ray and an edge
!    intersect and return the intersection point if they do.
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
!    Input, integer ( kind = 4 ) MODE, 0 for two edges, 1 (or nonzero) for a ray and an edge.
!
!    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, XW1, YW1, XW2, YW2,
!    the vertex coordinates; an edge (ray) is from (XV1,YV1) to (thru)
!    (XV2,YV2); an edge joins vertices (XW1,YW1) and (XW2,YW2).
!
!    Output, real ( kind = 8 ) XU, YU, the coordinates of the point of
!    intersection iff INTSCT is TRUE.
!
!    Output, logical INTSCT, TRUE if the edges/ray are nondegenerate, not
!    parallel, and intersect, FALSE otherwise.
!
  implicit none

  real    ( kind = 8 ) denom
  real    ( kind = 8 ) dxv
  real    ( kind = 8 ) dxw
  real    ( kind = 8 ) dyv
  real    ( kind = 8 ) dyw
  logical              intsct
  integer ( kind = 4 ) mode
  real    ( kind = 8 ) t
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
  intsct = .false.
  dxv = xv2 - xv1
  dyv = yv2 - yv1
  dxw = xw2 - xw1
  dyw = yw2 - yw1
  tolabs = tol * max ( abs ( dxv ), abs ( dyv ), abs ( dxw ), abs ( dyw ) )
  denom = dyv * dxw - dxv * dyw

  if ( abs ( denom ) <= tolabs) then
    return
  end if

  t = ( dyv * ( xv1 - xw1 ) - dxv * ( yv1 - yw1 ) ) / denom

  if ( t < -tol .or. 1.0D+00 + tol < t ) then
    return
  end if

  xu = xw1 + t * dxw
  yu = yw1 + t * dyw

  if ( abs ( dyv ) <= abs ( dxv ) ) then
    t = (xu - xv1)/dxv
  else
    t = ( yu - yv1 ) / dyv
  end if

  if ( mode == 0 ) then
    if ( -tol <= t .and. t <= 1.0D+00 + tol ) then
      intsct = .true.
    end if
  else
    if ( -tol <= t ) then
      intsct = .true.
    end if
  end if

  return
end
