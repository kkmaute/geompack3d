subroutine hexagon_vertices_2d ( x, y )

!*****************************************************************************80
!
!! HEXAGON_VERTICES_2D returns the vertices of the unit hexagon in 2D.
!
!  Diagram:
!
!      120_____60
!        /     \
!    180/       \0
!       \       /
!        \_____/
!      240     300
!
!  Discussion:
!
!    The unit hexagon has maximum radius 1, and is the hull of the points
!
!      (   1,              0 ),
!      (   0.5,   sqrt (3)/2 ),
!      ( - 0.5,   sqrt (3)/2 ),
!      ( - 1,              0 ),
!      ( - 0.5, - sqrt (3)/2 ),
!      (   0.5, - sqrt (3)/2 ).
!
!  Modified:
!
!    21 September 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, real ( kind = 8 ) X(6), Y(6), the coordinates of the vertices.
!
  implicit none

  real    ( kind = 8 ), parameter :: a = 0.8660254037844386D+00
  real    ( kind = 8 ) x(6)
  real    ( kind = 8 ) y(6)

  x(1:6) = (/ 1.0D+00, 0.5D+00, -0.5D+00, -1.0D+00, -0.5D+00,  0.5D+00 /)
  y(1:6) = (/ 0.0D+00, a,        a,        0.0D+00, -a,       -a /)

  return
end
