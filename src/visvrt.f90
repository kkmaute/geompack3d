subroutine visvrt ( angspc, xeye, yeye, nvis, xc, yc, ivis, maxn, nvsvrt, &
  theta )

!*****************************************************************************80
!
!! VISVRT determines a list of visible vertices.
!
!  Discussion:
!
!    This routine determines a list of visible vertices, ordered by
!    increasing "polar angle", on the boundary of the visibility
!    polygon from boundary eyepoint (XEYE,YEYE). This list
!    includes the vertices of visibility polygon such that a
!    line segment from (XEYE,YEYE) to vertex lies in interior
!    of polygon, as well as extra points on edges which subtend
!    an angle greater than or equal to 2 * ANGSPC at (XEYE,YEYE).
!
!    These extra points are at an equal angular spacing which is at
!    least ANGSPC and less than 2 * ANGSPC. The successor and predecessor
!    of eyepoint are included in list.
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
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter in
!    radians which controls how many extra points become visible vertices.
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of boundary eyepoint.
!
!    Input, integer ( kind = 4 ) NVIS, (number of vertices of visibility polygon) - 2.
!
!    Input/output, real ( kind = 8 ) XC(0:NVIS),YC(0:NVIS).  On input, the
!    coordinates of the vertices of visibility polygon in counterclockwise
!    order; (XC(0),YC(0)) and (XC(NVIS),YC(NVIS)) are the successor and
!    predecessor vertices of eyepoint in visibility polygon; at most 2
!    consecutive vertices have same polar angle with respect to
!    eyepoint.
!    On output, XC(0:NVSVRT),YC(0:NVSVRT) contain coordinates of
!    visible vertices which overwrite the input coordinates
!
!    Input/output, integer ( kind = 4 ) IVIS(0:NVIS).  On input, contains information
!    about the vertices of XC, YC arrays with respect to the original
!    polygon from which visibility polygon is computed; if 0 <= IVIS(I)
!    then (XC(I),YC(I)) has index I in original polygon;
!    if IVIS(I) < 0 then (XC(I),YC(I)) is on the edge
!    ending at vertex of index -IVIS(I) in original polygon;
!    indexing starts at 0 from successor of eyepoint.
!    On output, IVIS(0:NVSVRT) contains information about the output
!    vertices of XC, YC arrays as described above for input
!
!    Input, integer ( kind = 4 ) MAXN, the upper bound on NVSVRT; should be at least
!    NVIS + INT(PHI/ANGSPC) where PHI is the interior
!    angle at (XEYE,YEYE).
!
!    Output, integer ( kind = 4 ) NVSVRT, (number of visible vertices) - 1.
!
!    Output, real ( kind = 8 ) THETA(0:NVSVRT), the polar angles of visible
!    vertices with respect to (XEYE, YEYE) at origin and (XC(0),YC(0))
!    on positive x-axis.
!
  implicit none

  integer ( kind = 4 ) maxn

  real    ( kind = 8 ) alpha
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) ang1
  real    ( kind = 8 ) ang2
  real    ( kind = 8 ) angdif
  real    ( kind = 8 ) angle
  real    ( kind = 8 ) angsp2
  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) cosang
  integer ( kind = 4 ) cur
  real    ( kind = 8 ) diff
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) ivis(0:maxn)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) n
  real    ( kind = 8 ) numer
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvsvrt
  real    ( kind = 8 ) r
  real    ( kind = 8 ) sinang
  real    ( kind = 8 ) theta(0:maxn)
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:maxn)
  real    ( kind = 8 ) xeye
  real    ( kind = 8 ) yc(0:maxn)
  real    ( kind = 8 ) yeye
!
!  Shift input vertices right, and possibly remove first and last
!  vertices due to collinearity with eyepoint.
!
  tol = 100.0D+00 * epsilon ( tol )
  angsp2 = 2.0D+00 * angspc
  cur = maxn + 1
  n = maxn

  do i = nvis, 0, -1
    cur = cur - 1
    xc(cur) = xc(i)
    yc(cur) = yc(i)
    ivis(cur) = ivis(i)
  end do

  lr = lrline(xc(cur+1),yc(cur+1),xeye,yeye,xc(cur),yc(cur),0.0D+00)

  if ( 0 <= lr ) then
    cur = cur + 1
    xc(0) = xc(cur)
    yc(0) = yc(cur)
    ivis(0) = ivis(cur)
  end if

  lr = lrline(xc(n-1),yc(n-1),xeye,yeye,xc(n),yc(n),0.0D+00)

  if ( lr <= 0 ) then
    n = n - 1
  end if

  alpha = atan2 ( yc(0)-yeye, xc(0)-xeye )
  ang2 = 0.0D+00
  theta(0) = 0.0D+00
  top = 0
  cur = cur + 1
!
!  Process edge from vertices of indices CUR-1, CUR.
!
20 continue

  ang1 = ang2
  ang2 = angle(xc(cur),yc(cur),xeye,yeye,xc(0),yc(0))
  angdif = ang2 - ang1

  if ( angdif <= tol ) then

    diff = ((xc(cur) - xeye)**2 + (yc(cur) - yeye)**2) - &
      ((xc(cur-1) - xeye)**2 + (yc(cur-1) - yeye)**2)

    if ( diff < 0.0D+00 ) then
      xc(top) = xc(cur)
      yc(top) = yc(cur)
      ivis(top) = ivis(cur)
      theta(top) = ang2
    end if

  else

    if ( angsp2 <= angdif ) then

      k = int ( angdif / angspc )
      ind = -abs(ivis(cur))
      angdif = angdif / real ( k, kind = 8 )
      dx = xc(cur) - xc(cur-1)
      dy = yc(cur) - yc(cur-1)
      numer = (xc(cur) - xeye)*dy - (yc(cur) - yeye)*dx

      do i = 1, k-1
        top = top + 1
        theta(top) = ang1 + real ( i, kind = 8 ) * angdif
        ang = theta(top) + alpha
        cosang = cos(ang)
        sinang = sin(ang)
        r = numer/(dy*cosang - dx*sinang)
        xc(top) = r*cosang + xeye
        yc(top) = r*sinang + yeye
        ivis(top) = ind
      end do

    end if

    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = ivis(cur)
    theta(top) = ang2

  end if

  cur = cur + 1

  if ( cur <= n ) then
    go to 20
  end if

  nvsvrt = top

  return
end
