subroutine walkt2 ( x, y, ntri, vcl, til, tnbr, itri, iedg, ierr )

!*****************************************************************************80
!
!! WALKT2 walks through a 2D triangulation searching for a point.
!
!  Discussion:
!
!    This routine walks through neighboring triangles of a 2D (Delaunay)
!    triangulation until a triangle is found containing point (X,Y)
!    or (X,Y) is found to be outside the convex hull.  The search is
!    guaranteed to terminate for a Delaunay triangulation, else a
!    cycle may occur.
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
!    Input, real ( kind = 8 ) X, Y, the 2D point.
!
!    Input, integer ( kind = 4 ) NTRI, the number of triangles in triangulation;
!    used to detect cycle.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list.
!
!    Input/output, integer ( kind = 4 ) ITRI.  On input, the index of triangle to begin
!    search at.  On output, the index of triangle that search ends at.
!
!    Output, integer ( kind = 4 ) IEDG, 0 if ( X,Y) is in the interior of triangle ITRI;
!    I = 1, 2, or 3 if ( X,Y) is on interior of edge I of ITRI;
!    I = 4, 5, or 6 if ( X,Y) is (nearly) vertex I-3 of ITRI;
!    I = -1, -2, or -3 if ( X,Y) is outside convex hull due
!    to walking past edge -I of triangle ITRI.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) a
  real    ( kind = 8 ) alpha
  integer ( kind = 4 ) b
  real    ( kind = 8 ) beta
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cnt
  real    ( kind = 8 ) det
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dxa
  real    ( kind = 8 ) dxb
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dya
  real    ( kind = 8 ) dyb
  real    ( kind = 8 ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) itri
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Use barycentric coordinates to determine where (X,Y) is located.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  cnt = 0
  iedg = 0

10 continue

  cnt = cnt + 1

  if ( ntri < cnt ) then
    ierr = 226
    return
  end if

  a = til(1,itri)
  b = til(2,itri)
  c = til(3,itri)
  dxa = vcl(1,a) - vcl(1,c)
  dya = vcl(2,a) - vcl(2,c)
  dxb = vcl(1,b) - vcl(1,c)
  dyb = vcl(2,b) - vcl(2,c)
  dx = x - vcl(1,c)
  dy = y - vcl(2,c)
  det = dxa*dyb - dya*dxb
  alpha = (dx*dyb - dy*dxb) / det
  beta = (dxa*dy - dya*dx) / det
  gamma = 1.0D+00 - alpha - beta

  if ( tol < alpha .and. tol < beta .and. tol < gamma ) then

    return

  else if ( alpha < -tol ) then

    i = tnbr(2,itri)
    if ( i <= 0 ) then
      iedg = -2
      return
    end if

  else if ( beta < -tol ) then

    i = tnbr(3,itri)
    if ( i <= 0 ) then
      iedg = -3
      return
    end if

  else if ( gamma < -tol ) then

    i = tnbr(1,itri)
    if ( i <= 0 ) then
      iedg = -1
      return
    end if

  else if ( alpha <= tol ) then

    if ( beta <= tol ) then
      iedg = 6
    else if ( gamma <= tol ) then
      iedg = 5
    else
      iedg = 2
    end if
    return

  else if ( beta <= tol ) then

    if ( gamma <= tol ) then
      iedg = 4
    else
      iedg = 3
    end if
    return

  else

    iedg = 1
    return

  end if

  itri = i
  go to 10

end
