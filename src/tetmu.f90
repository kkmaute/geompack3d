function tetmu ( crit, a, b, c, d, s )

!*****************************************************************************80
!
!! TETMU computes a tetrahedron shape measure.
!
!  Discussion:
!
!    This routine computes a tetrahedron shape measure, scaled so that
!    the largest possible value is 1.
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
!    Input, integer ( kind = 4 ) CRIT, the criterion number:
!    1 for sin(min solid angle / 2),
!    2 for radius ratio,
!    3 for mean ratio.
!
!    Input, real ( kind = 8 ) A(1:3), B(1:3), C(1:3), D(1:3), the vertices
!    of the tetrahedron.
!
!    Workspace, real ( kind = 8 ) S(1:4), needed for SANGMN routine
!
!    Output, real ( kind = 8 ) TETMU, depending on CRIT, contains:
!    1, SANGMN(ABCD) * SF,
!    2, RADRTH(ABCD),
!    3, EMNRTH(ABCD).
!
  implicit none

  real    ( kind = 8 ) a(3)
  real    ( kind = 8 ) b(3)
  real    ( kind = 8 ) c(3)
  integer ( kind = 4 ) crit
  real    ( kind = 8 ) d(3)
  real    ( kind = 8 ) emnrth
  real    ( kind = 8 ) radrth
  real    ( kind = 8 ) s(4)
  real    ( kind = 8 ) sangmn
  real    ( kind = 8 ), parameter :: sf = 3.674234614174767D+00
  real    ( kind = 8 ) tetmu

  if ( crit == 1 ) then
    tetmu = sangmn ( a, b, c, d, s ) * sf
  else if ( crit == 2 ) then
    tetmu = radrth ( a, b, c, d )
  else if ( crit == 3 ) then
    tetmu = emnrth ( a, b, c, d )
  else
    tetmu = 1.0D+00
  end if

  return
end
