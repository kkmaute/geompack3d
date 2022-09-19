function opside ( a, b, c, d, e )

!*****************************************************************************80
!
!! OPSIDE tests if points are on opposite sides of a triangular face.
!
!  Discussion:
!
!    This routine tests if points D, E are on opposite sides of triangular
!    face with vertices A, B, C.
!
!  Modified:
!
!    31 August 2005
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), E(1:3),
!    five 3D points.
!
!    Output, integer ( kind = 4 ) OPSIDE, the result of the test:
!    +1 if D, E on opposite sides;
!    -1 if on same side;
!     2 if D is coplanar with face ABC (ABCD is a degenerate tetrahedron);
!     0 if E is coplanar with face ABC
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) c(3)
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) ddp
  real    ( kind = 8 ) dmax
  real    ( kind = 8 ) e(3)
  real    ( kind = 8 ) edp
  real    ( kind = 8 ) emax
  real    ( kind = 8 ) nrml1
  real    ( kind = 8 ) nrml2
  real    ( kind = 8 ) nrml3
  integer ( kind = 4 ) opside
  real    ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)

  emax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ) )

  dmax = max ( emax, abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  nrml1 = ab(2) * ac(3) - ab(3) * ac(2)
  nrml2 = ab(3) * ac(1) - ab(1) * ac(3)
  nrml3 = ab(1) * ac(2) - ab(2) * ac(1)

  ddp = ( d(1) - a(1) ) * nrml1 &
      + ( d(2) - a(2) ) * nrml2 &
      + ( d(3) - a(3) ) * nrml3

  if ( abs ( ddp ) <= tol * dmax ) then
    opside = 2
    return
  end if

  emax = max ( emax, abs ( e(1) ), abs ( e(2) ), abs ( e(3) ) )

  edp = ( e(1) - a(1) ) * nrml1 &
      + ( e(2) - a(2) ) * nrml2 &
      + ( e(3) - a(3) ) * nrml3

  if ( abs ( edp ) <= tol * emax ) then
    opside = 0
  else if ( ddp * edp < 0.0D+00 ) then
    opside = 1
  else
    opside = -1
  end if

  return
end
