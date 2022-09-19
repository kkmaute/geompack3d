subroutine shrnk2 ( nvrt, xc, yc, sdist, nshr, xs, ys, iedge, ierr )

!*****************************************************************************80
!
!! SHRNK2 shrinks a convex polygon.
!
!  Discussion:
!
!    This routine shrinks a convex polygon, with vertices given in
!    counterclockwise order and with all interior angles < PI, by distance
!    SDIST(I) for Ith edge, I = 0,...,NVRT-1.
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
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates
!    in counterclockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) SDIST(0:NVRT-1), the nonnegative shrink
!    distances for edges.
!
!    Output, integer ( kind = 4 ) NSHR, the number of vertices on boundary of shrunken
!    polygon; 0 if shrunken polygon is empty else 3 <= NSHR <= NVRT.
!
!    Output, real ( kind = 8 ) XS(0:NSHR), YS(0:NSHR), the coordinates
!    of shrunken polygon in counterclockwise order if NSHR is greater than 0;
!    (XS(0),YS(0)) = (XS(NSHR),YS(NSHR)).
!
!    Output, integer ( kind = 4 ) IEDGE(0:NSHR), the indices of edges of original polygon
!    in range 0 to NVRT-1 corresponding to each shrunken polygon edge.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) alpha
  logical              first
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge(0:nvrt)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nshr
  logical              parall
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  real    ( kind = 8 ) piptol
  real    ( kind = 8 ) sdist(0:nvrt-1)
  real    ( kind = 8 ) theta
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xs(0:nvrt)
  real    ( kind = 8 ) y
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) ys(0:nvrt)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pi2 = 2.0D+00 * pi
  piptol = pi + tol
  alpha = atan2 ( yc(1)-yc(0), xc(1)-xc(0) )
  i = 1
  call xline(xc(0),yc(0),xc(1),yc(1),xc(1),yc(1),xc(2),yc(2), &
    sdist(0),sdist(1),xs(1),ys(1),parall)

  if ( parall ) then
    if ( abs(sdist(0)-sdist(1)) <= max ( tol, 1.0D-05 ) * sdist(0) ) then
      i = 2
      call xline(xc(0),yc(0),xc(1),yc(1),xc(2),yc(2),xc(3),yc(3), &
        sdist(0),sdist(2),xs(1),ys(1),parall)
    end if
  end if

  if ( parall ) then
    ierr = 202
    nshr = 0
    return
  end if

  iedge(0) = 0
  iedge(1) = i
  i = i + 1
  j = 0
  nshr = 1
  first = .true.
!
!  First while loop processes edges subtending angle <= PI
!  with respect to first edge.
!
  do

    theta = atan2 ( yc(i+1)-yc(i), xc(i+1)-xc(i) ) - alpha

    if ( theta < 0.0D+00 ) then
      theta = theta + pi2
    end if

    if ( piptol < theta ) then
      call xline ( xc(0), yc(0), xc(1), yc(1), xc(i), yc(i), xc(i+1), &
        yc(i+1), 0.0D+00, 0.0D+00, x, y, parall )
      if ( .not. parall ) then
        go to 40
      end if
    end if

20  continue

    lr = lrline(xs(nshr),ys(nshr),xc(i),yc(i),xc(i+1),yc(i+1),sdist(i))

    if ( lr < 0 ) then
      go to 30
    end if

    nshr = nshr - 1

    if ( 1 <= nshr ) then
      go to 20
    end if

30  continue

    if ( nshr < 1 .and. abs ( theta - pi ) <= tol ) then
      nshr = 0
      return
    end if

    k = iedge(nshr)
    nshr = nshr + 1
    call xline(xc(k),yc(k),xc(k+1),yc(k+1),xc(i),yc(i),xc(i+1), &
      yc(i+1),sdist(k),sdist(i),xs(nshr),ys(nshr),parall)

    if ( parall ) then
      nshr = 0
      return
    end if

    iedge(nshr) = i
    i = i + 1

  end do
!
!  Second while loop processes remaining edges.
!
40 continue

  if ( first ) then
    first = .false.
    go to 50
  end if

  lr = lrline(xs(j),ys(j),xc(i),yc(i),xc(i+1),yc(i+1),sdist(i))

  if ( lr <= 0) go to 70

50 continue

  if ( nshr <= j ) then
    nshr = 0
    return
  end if

  lr = lrline(xs(nshr),ys(nshr),xc(i),yc(i),xc(i+1), &
    yc(i+1),sdist(i))

  if ( 0 <= lr ) then
    nshr = nshr - 1
    go to 50
  end if

  k = iedge(nshr)
  nshr = nshr + 1
  call xline(xc(k),yc(k),xc(k+1),yc(k+1),xc(i),yc(i),xc(i+1), &
    yc(i+1),sdist(k),sdist(i),xs(nshr),ys(nshr),parall)

  if ( parall ) then
    ierr = 202
    nshr = 0
    return
  end if

  iedge(nshr) = i

60 continue

  lr = lrline(xs(j+1),ys(j+1),xc(i),yc(i),xc(i+1),yc(i+1), sdist(i))

  if ( 0 <= lr ) then
    j = j + 1
    go to 60
  end if

  k = iedge(j)
  call xline(xc(k),yc(k),xc(k+1),yc(k+1),xc(i),yc(i),xc(i+1), &
    yc(i+1),sdist(k),sdist(i),xs(j),ys(j),parall)

  if ( parall ) then
    ierr = 202
    nshr = 0
    return
  end if

  xs(nshr+1) = xs(j)
  ys(nshr+1) = ys(j)
  iedge(nshr+1) = iedge(j)

70 continue

  i = i + 1
  if ( i < nvrt) go to 40

  if ( 0 < j ) then
    do i = 0, nshr+1-j
      xs(i) = xs(i+j)
      ys(i) = ys(i+j)
      iedge(i) = iedge(i+j)
    end do
  end if

  nshr = nshr + 1 - j

  return
end
