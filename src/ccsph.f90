subroutine ccsph ( intest, a, b, c, d, e, center, radsq, in )

!*****************************************************************************80
!
!! CCSPH finds the circumsphere through the vertices of a tetrahedron.
!
!  Discussion:
!
!    This routine finds the center and the square of the radius of
!    the circumsphere through four vertices of a tetrahedron, and
!    possibly determines whether a fifth 3D point is inside the sphere.
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
!    Input, logical INTEST, is TRUE, iff test for fifth point in sphere
!    to be made.
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), vertices
!    of tetrahedron.
!
!    Input, real ( kind = 8 ) E(1:3), a fifth point; referenced iff
!    INTEST is TRUE.
!
!    Output, real ( kind = 8 ) CENTER(1:3), center of sphere; undefined
!    if A,B,C,D coplanar.
!
!    Output, real ( kind = 8 ) RADSQ, the square of radius of sphere;
!    -1 if A,B,C,D coplanar.
!
!    Output, integer ( kind = 4 ) IN, contains following value if INTEST is .TRUE.:
!     2 if A,B,C,D coplanar;
!     1 if E inside sphere;
!     0 if E on sphere;
!    -1 if E outside sphere
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) c(3)
  real    ( kind = 8 ) center(3)
  real    ( kind = 8 ) cmax
  real    ( kind = 8 ) cp1
  real    ( kind = 8 ) cp2
  real    ( kind = 8 ) cp3
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) da(3)
  real    ( kind = 8 ) db(3)
  real    ( kind = 8 ) dc(3)
  real    ( kind = 8 ) det
  real    ( kind = 8 ) dsq
  real    ( kind = 8 ) e(3)
  integer ( kind = 4 ) in
  logical              intest
  real    ( kind = 8 ) radsq
  real    ( kind = 8 ) rhs(3)
  real    ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )

  da(1:3) = a(1:3) - d(1:3)
  db(1:3) = b(1:3) - d(1:3)
  dc(1:3) = c(1:3) - d(1:3)

  rhs(1) = 0.5D+00 * sum ( da(1:3)**2 )
  rhs(2) = 0.5D+00 * sum ( db(1:3)**2 )
  rhs(3) = 0.5D+00 * sum ( dc(1:3)**2 )

  cmax = max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ), &
    abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  cp1 = db(2) * dc(3) - dc(2) * db(3)
  cp2 = dc(2) * da(3) - da(2) * dc(3)
  cp3 = da(2) * db(3) - db(2) * da(3)

  det = da(1) * cp1 + db(1) * cp2 + dc(1) * cp3

  if ( abs ( det ) <= 0.01D+00 * tol * cmax ) then
    radsq = -1.0D+00
    in = 2
    return
  end if

  center(1) = ( rhs(1) * cp1 + rhs(2) * cp2 + rhs(3) * cp3 ) / det

  cp1 = db(1) * rhs(3) - dc(1) * rhs(2)
  cp2 = dc(1) * rhs(1) - da(1) * rhs(3)
  cp3 = da(1) * rhs(2) - db(1) * rhs(1)

  center(2) =  ( da(3) * cp1 + db(3) * cp2 + dc(3) * cp3 ) / det
  center(3) = -( da(2) * cp1 + db(2) * cp2 + dc(2) * cp3 ) / det

  radsq = sum ( center(1:3)**2 )

  center(1:3) = center(1:3) + d(1:3)

  if ( intest ) then

    dsq = sum ( ( e(1:3) - center(1:3) )**2 )

    if ( ( 1.0D+00 + tol ) * radsq < dsq ) then
      in = -1
    else if ( dsq < ( 1.0D+00 - tol ) * radsq ) then
      in = 1
    else
      in = 0
    end if

  end if

  return
end
