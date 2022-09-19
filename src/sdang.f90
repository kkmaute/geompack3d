subroutine sdang ( a, b, c, d, sang, dang )

!*****************************************************************************80
!
!! SDANG computes the solid and dihedral angles of a tetrahedron.
!
!  Discussion:
!
!    This routine computes the 4 solid angles (+ PI) and 6 dihedral angles of
!    a tetrahedron in radians.  The solid angle at vertex A is
!    (dihedral angle at edge AB) + (dihedral angle at edge AC)
!    + (dihedral angle at edge AD) - PI
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), the vertices
!    of tetrahedron.
!
!    Output, real ( kind = 8 ) SANG(1:4), the solid angles (+ PI) at
!    A,B,C,D, respectively.
!
!    Output, real ( kind = 8 ) DANG(1:6), the dihedral angles at
!    AB,AC,AD,BC,BD,CD, respectively.
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ad(3)
  real    ( kind = 8 ) angle3
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) bc(3)
  real    ( kind = 8 ) bd(3)
  real    ( kind = 8 ) c(3)
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) dang(6)
  real    ( kind = 8 ) nrmabc(3)
  real    ( kind = 8 ) nrmabd(3)
  real    ( kind = 8 ) nrmacd(3)
  real    ( kind = 8 ) nrmbcd(3)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) rtolsq
  real    ( kind = 8 ) sang(4)
  real    ( kind = 8 ) tol
!
!  Compute outward (or inward) normals at 4 faces to get angles.
!
  tol = 100.0D+00 * epsilon ( tol )

  ab(1:3) = b(1:3) - a(1:3)
  ac(1:3) = c(1:3) - a(1:3)
  ad(1:3) = d(1:3) - a(1:3)
  bc(1:3) = c(1:3) - b(1:3)
  bd(1:3) = d(1:3) - b(1:3)

  rtolsq = tol * max ( &
    abs ( a(1) ), abs ( a(2) ), abs ( a(3) ), &
    abs ( b(1) ), abs ( b(2) ), abs ( b(3) ), &
    abs ( c(1) ), abs ( c(2) ), abs ( c(3) ), &
    abs ( d(1) ), abs ( d(2) ), abs ( d(3) ) )

  rtolsq = rtolsq**2

  nrmabc(1) = ac(2) * ab(3) - ac(3) * ab(2)
  nrmabc(2) = ac(3) * ab(1) - ac(1) * ab(3)
  nrmabc(3) = ac(1) * ab(2) - ac(2) * ab(1)
  nrmabd(1) = ab(2) * ad(3) - ab(3) * ad(2)
  nrmabd(2) = ab(3) * ad(1) - ab(1) * ad(3)
  nrmabd(3) = ab(1) * ad(2) - ab(2) * ad(1)
  nrmacd(1) = ad(2) * ac(3) - ad(3) * ac(2)
  nrmacd(2) = ad(3) * ac(1) - ad(1) * ac(3)
  nrmacd(3) = ad(1) * ac(2) - ad(2) * ac(1)
  nrmbcd(1) = bc(2) * bd(3) - bc(3) * bd(2)
  nrmbcd(2) = bc(3) * bd(1) - bc(1) * bd(3)
  nrmbcd(3) = bc(1) * bd(2) - bc(2) * bd(1)

  dang(1) = pi - angle3 ( nrmabc, nrmabd, rtolsq )
  dang(2) = pi - angle3 ( nrmabc, nrmacd, rtolsq )
  dang(3) = pi - angle3 ( nrmabd, nrmacd, rtolsq )
  dang(4) = pi - angle3 ( nrmabc, nrmbcd, rtolsq )
  dang(5) = pi - angle3 ( nrmabd, nrmbcd, rtolsq )
  dang(6) = pi - angle3 ( nrmacd, nrmbcd, rtolsq )

  sang(1) = dang(1) + dang(2) + dang(3)
  sang(2) = dang(1) + dang(4) + dang(5)
  sang(3) = dang(2) + dang(4) + dang(6)
  sang(4) = dang(3) + dang(5) + dang(6)

  return
end
