subroutine sepshp ( angtol, nvrt, xc, yc, indpvl, iang, i1, i2, wk, ierr )

!*****************************************************************************80
!
!! SEPSHP determines a separator that splits a convex polygon.
!
!  Discussion:
!
!    This routine determines a separator to split convex polygon into two
!    parts based on shape (diameter) of polygon.
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
!    Input, real ( kind = 8 ) ANGTOL, angle tolerance parameter (in radians).
!
!    Input, integer ( kind = 4 ) NVRT, number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), coordinates of polygon
!    vertices in counterclockwise order, translated so that centroid is at
!    origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, integer ( kind = 4 ) INDPVL(0:NVRT), indices in PVL of vertices;
!    INDPVL(I) = -K if ( XC(I),YC(I)) is extra vertex inserted on edge from
!    K to PVL(SUCC,K).
!
!    Input, real ( kind = 8 ) IANG(1:*), interior angle array
!
!    Output, integer ( kind = 4 ) I1, I2, indices in range 0 to NVRT-1 of best separator
!    according to shape and max-min angle criterion; I1 = -1
!    if no satisfactory separator is found.
!
!    Workspace, real ( kind = 8 ) WK(1:2*NVRT).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indpvl(0:nvrt)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pimtol
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v(2)
  integer ( kind = 4 ) w(2)
  real    ( kind = 8 ) wk(2*nvrt)
  real    ( kind = 8 ) xa
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) ya
  real    ( kind = 8 ) yc(0:nvrt)
!
!  Determine diameter of polygon.  Possible separators endpoints (two
!  on each side of polygon) are nearest to perpendicular bisector of
!  diameter.  (XA,YA) and (XA+DX,YA+DY) are on bisector.  Distance to
!  bisector is proportional to two times triangle area.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pimtol = pi - tol
  n = 0

  do i = 0, nvrt-1
    k = indpvl(i)
    if ( 0 < k ) then
      if ( iang(k) < pimtol ) then
        n = n + 1
        wk(n) = xc(i)
        wk(n+nvrt) = yc(i)
      end if
    end if
  end do

  call diam2 ( n, wk, wk(nvrt+1), i1, i2, dist, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  if ( i2 < i1 ) then
    i = i1
    i1 = i2
    i2 = i
  end if

  dx = wk(i2+nvrt) - wk(i1+nvrt)
  dy = wk(i1) - wk(i2)
  xa = 0.5D+00 * ( wk(i1) + wk(i2) - dx )
  ya = 0.5D+00 * ( wk(i1+nvrt) + wk(i2+nvrt) - dy )

  i = i1 - 1

  do

    if ( xc(i) == wk(i1) .and. yc(i) == wk(i1+nvrt) ) then
      i1 = i
      exit
    end if

    i = i + 1

  end do

  i = max ( i2-1, i1+1 )

  do

    if ( xc(i) == wk(i2) .and. yc(i) == wk(i2+nvrt) ) then
      i2 = i
      exit
    end if

    i = i + 1

  end do

  i = i1 + 1

  do

    dist = dx * ( yc(i) - ya ) - dy * ( xc(i) - xa )

    if ( 0.0D+00 <= dist ) then
      v(1) = i - 1
      v(2) = i
      exit
    end if

    i = i + 1

  end do

  i = i2 + 1

  do

    if ( nvrt <= i ) then
      i = 0
    end if

    dist = dx * ( yc(i) - ya) - dy * ( xc(i) - xa)

    if ( dist <= 0.0D+00 ) then
      w(1) = i - 1
      w(2) = i
      if ( i <= 0 ) then
        w(1) = nvrt - 1
      end if
      exit
    end if

    i = i + 1

  end do

  call mmasep ( angtol, xc, yc, indpvl, iang, v, w, i1, i2 )

  return
end
