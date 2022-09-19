subroutine vornbr ( xeye, yeye, nvrt, xc, yc, nvor, ivor, xvor, yvor, ierr )

!*****************************************************************************80
!
!! VORNBR determines the Voronoi neighbors of a point.
!
!  Discussion:
!
!    This routine determines the Voronoi neighbors of (XEYE,YEYE) from a
!    list of vertices which are in increasing "polar angle" order.
!
!    The Voronoi neighbors are a sublist of this list. The
!    Voronoi polygon is restricted to the sector formed from the
!    the edges joining (XEYE,YEYE) to the first and last vertices
!    of this list. Each Voronoi neighbor corresponds to an edge
!    of the Voronoi polygon.
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
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of eyepoint.
!
!    Input, integer ( kind = 4 ) NVRT, (number of vertices in list) - 1.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates
!    from which Voronoi neighbors are determined; (XC(0),YC(0)),...,
!    (XC(NVRT),YC(NVRT)) are in increasing angular
!    displacement order with respect to (XEYE,YEYE).
!
!    Output, integer ( kind = 4 ) NVOR, (number of Voronoi neighbors) - 1 [<= NVRT].
!
!    Output, integer ( kind = 4 ) IVOR(0:NVOR), the indices of Voronoi neighbors in XC, YC
!    arrays; 0 <= IVOR(0) < ... < IVOR(NVOR) <= NVRT.
!
!    Workspace, real ( kind = 8 ) XVOR(0:NVRT), YVOR(0:NVRT), arrays for
!    storing the vertex coordinates of the Voronoi polygon.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) a11
  real    ( kind = 8 ) a12
  real    ( kind = 8 ) a21
  real    ( kind = 8 ) a22
  real    ( kind = 8 ) b1
  real    ( kind = 8 ) b2
  real    ( kind = 8 ) det
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) im
  integer ( kind = 4 ) ivor(0:nvrt)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvor
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xeye
  real    ( kind = 8 ) xi
  real    ( kind = 8 ) xvor(0:nvrt)
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) yeye
  real    ( kind = 8 ) yi
  real    ( kind = 8 ) yvor(0:nvrt)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  k = 1
  m = 0
  ivor(0) = 0
  xvor(0) = ( xeye + xc(0) ) * 0.5D+00
  yvor(0) = ( yeye + yc(0) ) * 0.5D+00
!
!  Beginning of main loop
!
  do while ( k <= nvrt)
!
!  Determine the intersection of the perpendicular bisectors
!  of edges from (XEYE,YEYE) to (XC(K),YC(K)) and from
!  (XEYE,YEYE) to (XC(IM),YC(IM)).
!
    im = ivor(m)
    a11 = xc(k) - xeye
    a12 = yc(k) - yeye
    a21 = xc(im) - xeye
    a22 = yc(im) - yeye
    tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
    det = a11 * a22 - a21 * a12

    if ( abs(det) <= tolabs ) then
      ierr = 212
      return
    end if

    b1 = ( a11**2 + a12**2 ) * 0.5D+00
    b2 = ( a21**2 + a22**2 ) * 0.5D+00
    xi = ( b1 * a22 - b2 * a12 ) / det
    yi = ( b2 * a11 - b1 * a21 ) / det
!
!  Determine whether (XVOR(M+1),YVOR(M+1)) is to the left of or
!  on the directed line from (XEYE,YEYE) to (XVOR(M),YVOR(M)).
!
    xvor(m+1) = xi + xeye
    yvor(m+1) = yi + yeye
    lr = lrline(xvor(m+1),yvor(m+1),xeye,yeye,xvor(m),yvor(m), 0.0D+00)

    if ( lr <= 0 ) then
      m = m + 1
      ivor(m) = k
      k = k + 1
    else if ( 0 < m ) then
      m = m - 1
    else
!
!  Determine the intersection of edge from (XEYE,YEYE) to
!  (XC(0),YC(0)) and the perpendicular bisector of the edge
!  from (XEYE,YEYE) to (XC(K),YC(K)).
!
      a11 = xc(k) - xeye
      a12 = yc(k) - yeye
      a21 = yc(0) - yeye
      a22 = xeye - xc(0)
      tolabs = tol * max ( abs ( a11 ), abs ( a12 ), abs ( a21 ), abs ( a22 ) )
      det = a11 * a22 - a21 * a12

      if ( abs ( det ) <= tolabs ) then
        ierr = 212
        return
      end if

      b1 = ( a11**2 + a12**2 ) * 0.5D+00
      b2 = 0.0D+00
      xi = ( b1 * a22 - b2 * a12 ) / det
      yi = ( b2 * a11 - b1 * a21 ) / det
      xvor(m) = xi + xeye
      yvor(m) = yi + yeye
      ivor(m) = k
      k = k + 1

    end if

  end do
!
!  The following short loop determines which vertices at the end
!  of list are not Voronoi neighbors.
!
  do

    lr = lrline ( xvor(m), yvor(m), xeye, yeye, xc(nvrt), yc(nvrt), 0.0D+00 )

    if ( 0 <= lr ) then
      exit
    end if

    m = m - 1

    if ( m < 0 ) then
      exit
    end if

  end do

  nvor = m

  return
end
