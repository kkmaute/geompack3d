function mdf3 ( x, y, z, widp, nfcev, nedev, nvrev, &
  listev, infoev, ivrt, facval, edgval, vrtval, vcl )

!*****************************************************************************80
!
!! MDF3 evaluates a heuristic mesh distribution function in 3D.
!
!  Discussion:
!
!    This routine evaluates a heuristic mesh distribution function at (X,Y,Z).
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
!    Input, real ( kind = 8 ) X, Y, Z, the coordinates of 3D point.
!
!    Input, real ( kind = 8 ) WIDP, the width of polyhedron containing (X,Y,Z).
!
!    Input, integer ( kind = 4 ) NFCEV, NEDEV, NVREV, LISTEV(1:NFCEV+NEDEV+NVREV),
!    INFOEV(1:4,1:NFCEV+NEDEV), output from routine PRMDF3.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), real ( kind = 8 ) FACVAL(1:*), EDGVAL(1:*),
!    VRTVAL(1:*), arrays output from routine DSMDF3.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list
!
!    Output, real ( kind = 8 ) MDF3, the reciprocal of cube of length scale
!    at (X,Y,Z).
!
  implicit none

  real    ( kind = 8 ) cp(3)
  real    ( kind = 8 ) d
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) edgval(*)
  real    ( kind = 8 ) facval(*)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) infoev(4,*)
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev(*)
  real    ( kind = 8 ) mdf3
  integer ( kind = 4 ) nedev
  integer ( kind = 4 ) nfcev
  integer ( kind = 4 ) nvrev
  real    ( kind = 8 ) s
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) widp
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z

  s = widp
  k = 0

  do i = 1, nfcev
    k = k + 1
    d = abs ( infoev(1,k) * x + infoev(2,k) * y + infoev(3,k) * z - infoev(4,k))
    d = max ( 0.5D+00 * d, facval(listev(k)) )
    s = min ( s, d )
  end do

  do i = 1, nedev
    k = k + 1
    j = listev(k)
    l = ivrt(j)
    dir(1) = x - vcl(1,l)
    dir(2) = y - vcl(2,l)
    dir(3) = z - vcl(3,l)
    cp(1) = infoev(2,k)*dir(3) - infoev(3,k)*dir(2)
    cp(2) = infoev(3,k)*dir(1) - infoev(1,k)*dir(3)
    cp(3) = infoev(1,k)*dir(2) - infoev(2,k)*dir(1)
    d = sqrt(cp(1)**2 + cp(2)**2 + cp(3)**2)/infoev(4,k)
    d = max ( 0.5D+00 * d, edgval(j) )
    s = min ( s, d )
  end do

  do i = 1, nvrev
    k = k + 1
    j = listev(k)
    d = sqrt((vcl(1,j) - x)**2 + (vcl(2,j) - y)**2 + (vcl(3,j) - z)**2)
    d = max ( 0.5D+00 * d, vrtval(j) )
    s = min(s,d)
  end do

  mdf3 = 1.0D+00 / s**3

  return
end
