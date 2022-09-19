subroutine inttri ( nvrt, xc, yc, h, ibot, costh, sinth, ldv, nvc, ntri,  &
  maxvc, maxti, maxcw, vcl, til, ncw, cwalk, ierror )

!*****************************************************************************80
!
!! INTTRI generates triangles inside a convex polygon.
!
!  Discussion:
!
!    This routine generates triangles inside a convex polygon using
!    a quasi-uniform grid of spacing H.  It is assumed that the diameter
!    of the polygon is parallel to the Y axis.
!
!  Modified:
!
!    02 May 2001
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, real ( kind = 8 ) real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT),
!    the vertex coordinates in counterclockwise order;
!    (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) H, the spacing of mesh vertices in polygon.
!
!    Input, integer ( kind = 4 ) IBOT, the index of bottom vertex; diameter contains
!    vertices (XC(0),YC(0)) and (XC(IBOT),YC(IBOT)).
!
!    Input, real ( kind = 8 ) real ( kind = 8 ) COSTH, SINTH; COS(THETA),
!    SIN(THETA) where THETA in [-PI,PI] is rotation angle to get diameter
!    parallel to y-axis.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of coordinates or positions
!    used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NTRI, the number of triangles or positions
!    used in TIL.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXTI, the maximum size available for TIL array.
!
!    Input, integer ( kind = 4 ) MAXCW, the maximum size available for CWALK array;
!    assumed to be greater than or equal to 6*(1 + INT((YC(0) - YC(IBOT))/H)).
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Output, integer ( kind = 4 ) NCW, the number of mesh vertices in closed walk,
!    except NCW = 0 for 1 vertex.
!
!    Output, integer ( kind = 4 ) CWALK(0:NCW), indices in VCL of mesh vertices of closed
!    walk; CWALK(0) = CWALK(NCW)
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) maxcw
  integer ( kind = 4 ) maxti
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) a
  real    ( kind = 8 ) b
  real    ( kind = 8 ) costh
  integer ( kind = 4 ) cwalk(0:maxcw)
  real    ( kind = 8 ) cy
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) il
  integer ( kind = 4 ) im1l
  integer ( kind = 4 ) im1r
  integer ( kind = 4 ) ir
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) l0
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ncw
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  integer ( kind = 4 ) r
  integer ( kind = 4 ) r0
  integer ( kind = 4 ) r1
  integer ( kind = 4 ) rw
  real    ( kind = 8 ) sinth
  real    ( kind = 8 ) sy
  integer ( kind = 4 ) til(3,maxti)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(ldv,maxvc)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xj
  real    ( kind = 8 ) xk
  real    ( kind = 8 ) xl
  real    ( kind = 8 ) xm1l
  real    ( kind = 8 ) xm1r
  real    ( kind = 8 ) xr
  real    ( kind = 8 ) y
  real    ( kind = 8 ) yc(0:nvrt)

  ierror = 0

  tol = 100.0D+00 * epsilon ( tol )

  n = int ( ( yc(0) - yc(ibot) ) / h )
  y = yc(0) - 0.5D+00 * ( yc(0) - yc(ibot ) - real ( n, kind = 8 ) * h )
  l = 0
  r = nvrt

  do i = 0, n
!
!  Determine left and right x-coordinates of polygon for
!  scan line with y-coordinate Y, and generate mesh vertices.
!
    do while ( y < yc(l+1) )
      l = l + 1
    end do

    do while ( y < yc(r-1) )
      r = r - 1
    end do

    xl = xc(l) + ( xc(l+1) - xc(l) ) * ( y - yc(l) ) / ( yc(l+1) - yc(l) )
    xr = xc(r) + ( xc(r-1) - xc(r) ) * ( y - yc(r) ) / ( yc(r-1) - yc(r) )
    m = int ( ( xr - xl ) / h )
    x = xl + 0.5D+00 * ( xr - xl - real ( m, kind = 8 ) * h )

    if ( maxvc < nvc + m + 1 ) then
      ierror = 3
      return
    end if

    cy = costh * y
    sy = sinth * y
    il = nvc + 1
    xl = x

    do j = 0, m
      nvc = nvc + 1
      vcl(1,nvc) = costh * x + sy
      vcl(2,nvc) = cy - sinth * x
      x = x + h
    end do

    ir = nvc
    xr = x - h

    if ( n == 0 ) then

      ncw = 0
      cwalk(0) = nvc
      return

    else if ( i == 0 ) then

      lw = 0
      cwalk(lw) = il
      rw = maxcw + 1

      do j = il, ir
        rw = rw - 1
        cwalk(rw) = j
      end do

      go to 100

    end if
!
!  Generate triangles between scan lines Y+H and Y.
!
    a = max ( xl, xm1l )
    b = min ( xr, xm1r )

    if ( xm1l == a ) then
      l0 = im1l
      x = ( xm1l - xl ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l1 = il + j
    else
      l1 = il
      x = ( xl - xm1l ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      l0 = im1l + j
    end if

    if ( xm1r == b ) then
      r0 = im1r
      x = ( xr - xm1r ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r1 = ir - j
    else
      r1 = ir
      x = ( xm1r - xr ) / h
      j = int ( x + tol )
      if ( abs ( x - real ( j, kind = 8 ) ) <= tol ) then
        j = j - 1
      end if
      if ( j < 0 ) then
        j = 0
      end if
      r0 = im1r - j
    end if

    if ( l0 < r0 .or. l1 < r1 ) then

      j = l0
      k = l1
      xj = xm1l + real ( j-im1l, kind = 8 ) * h
      xk = xl + real ( k - il, kind = 8 ) * h

      do

        if ( k < r1 .and. ( xk <= xj .or. j == r0 ) ) then
          p = k
          k = k + 1
          xk = xk + h
        else
          p = j
          j = j + 1
          xj = xj + h
        end if

        ntri = ntri + 1

        if ( maxti < ntri ) then
          ierror = 9
          return
        end if

        til(1,ntri) = j
        til(2,ntri) = p
        til(3,ntri) = k

        if ( r0 <= j .and. r1 <= k ) then
          exit
        end if

      end do

    end if
!
!  Generate paths of closed walk between scan lines Y+H and Y.
!
    if ( xm1l < xl ) then
      do j = im1l+1, l0
        lw = lw + 1
        cwalk(lw) = j
      end do
      lw = lw + 1
      cwalk(lw) = il
    else
      do j = l1, il, -1
        lw = lw + 1
        cwalk(lw) = j
      end do
    end if

    if ( xr < xm1r ) then
      do j = im1r-1, r0, -1
        rw = rw - 1
        cwalk(rw) = j
      end do
      rw = rw - 1
      cwalk(rw) = ir
    else
      do j = r1, ir
        rw = rw - 1
        cwalk(rw) = j
      end do
    end if

100 continue

    y = y - h
    im1l = il
    im1r = ir
    xm1l = xl
    xm1r = xr

  end do
!
!  Add last path of left walk and shift indices of right walk.
!
  if ( m == 0 ) then
    rw = rw + 1
  else
    do j = il+1, ir-1
      lw = lw + 1
      cwalk(lw) = j
    end do
  end if

  if ( rw <= lw ) then
    ierror = 10
    return
  end if

  do j = rw, maxcw
    lw = lw + 1
    cwalk(lw) = cwalk(j)
  end do

  ncw = lw

  return
end
