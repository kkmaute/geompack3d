function volth ( a, b, c, d )

!*****************************************************************************80
!
!! VOLTH computes the volume of a tetrahedron.
!
!  Discussion:
!
!    This routine computes 6 times the volume of a tetrahedron.
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
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), vertices
!    of tetrahedron.
!
!    Output, real ( kind = 8 ) VOLTH, 6 * volume of tetrahedron.
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) c(3)
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) u(3)
  real    ( kind = 8 ) v(3)
  real    ( kind = 8 ) volth
  real    ( kind = 8 ) w(3)

  u(1:3) = b(1:3) - a(1:3)
  v(1:3) = c(1:3) - a(1:3)
  w(1:3) = d(1:3) - a(1:3)

  volth = abs( u(1)*(v(2)*w(3) - v(3)*w(2)) + u(2)*(v(3)*w(1) - &
    v(1)*w(3)) + u(3)*(v(1)*w(2) - v(2)*w(1)) )

  return
end
