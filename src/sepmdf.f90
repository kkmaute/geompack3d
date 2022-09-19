subroutine sepmdf ( angtol, nvrt, xc, yc, arpoly, mean, mdftr, indpvl, &
  iang, i1, i2 )

!*****************************************************************************80
!
!! SEPMDF determines a separator that splits a convex polygon.
!
!  Discussion:
!
!    This routine determines a separator to split a convex polygon into two
!    parts based on mesh distribution function.
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
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter
!    (in radians).
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the coordinates of
!    polygon vertices in counterclockwise order, translated so that
!    centroid is at origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) ARPOLY, the area of polygon.
!
!    Input, real ( kind = 8 ) MEAN, the mean mdf value in polygon.
!
!    Input, real ( kind = 8 ) MDFTR(0:NVRT-1), the mean mdf value in each
!    triangle of polygon; triangles are determined by polygon vertices
!    and centroid.
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), the indices in PVL of vertices;
!    INDPVL(I) = -K if (XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), the interior angle array
!
!    Output, integer ( kind = 4 ) I1, I2, the indices in range 0 to NVRT-1 of best separator
!    according to mdf and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) areatr
  real    ( kind = 8 ) arpoly
  integer ( kind = 4 ) hi
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) indpvl(0:nvrt)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mdftr(0:nvrt-1)
  real    ( kind = 8 ) mean
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) yc(0:nvrt)
!
!  Determine triangle with highest mean mesh density; then determine
!  triangles adjacent to this triangle with mesh density greater than
!  or equal to MEAN such that the area of these triangles is <= ARPOLY/2.
!  Note that twice the triangle area is computed.
!
  tol = 100.0D+00 * epsilon ( tol )
  hi = 0

  do i = 1, nvrt-1
    if ( mdftr(hi) < mdftr(i) ) then
      hi = i
    end if
  end do

  sum2 = xc(hi) * yc(hi+1) - xc(hi+1) * yc(hi)

  l = hi - 1
  if ( l < 0 ) then
    l = nvrt - 1
  end if

  m = hi + 1
  if ( nvrt <= m ) then
    m = 0
  end if

20 continue

  if ( mdftr(m) <= mdftr(l) ) then
    i = l
  else
    i = m
  end if

  if ( mdftr(i) < mean ) then
    go to 30
  end if

  areatr = xc(i) * yc(i+1) - xc(i+1) * yc(i)
  sum2 = sum2 + areatr
  if ( arpoly < sum2 ) go to 30

  if ( i == l ) then
    l = l - 1
    if ( l < 0 ) then
      l = nvrt - 1
    end if
  else
    m = m + 1
    if ( nvrt <= m ) then
      m = 0
    end if
  end if

  go to 20

30 continue

  l = l + 1
  if ( nvrt <= l ) then
    l = 0
  end if
!
!  Interchange role of L and M depending on angle determined by
!  (XC(M),YC(M)), (0,0), and (XC(L),YC(L)).
!  Possible separators are L,M; L,M+1; L+1,M; L+1,M+1.
!
  if ( pi < angle(xc(m),yc(m),0.0D+00,0.0D+00,xc(l),yc(l)) ) then
    i = l
    l = m
    m = i
  end if

  v(1) = l
  v(2) = l - 1
  if ( v(2) < 0 ) then
    v(2) = nvrt - 1
  end if
  w(1) = m
  w(2) = m + 1
  if ( nvrt <= w(2) ) then
    w(2) = 0
  end if

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
