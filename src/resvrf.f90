subroutine resvrf ( vr, aspc2d, atol2d, maxiw, maxwk, x, y, iang, link, w1, &
  w2, pt1, pt2, iwk, wk, ierr )

!*****************************************************************************80
!
!! RESVRF resolves a reflex vertex of a simple polygon.
!
!  Discussion:
!
!    This routine resolves a reflex vertex of simple polygon with one or two
!    separators, where polygon is face of polyhedron.
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
!    Input, integer ( kind = 4 ) VR, the index in X,Y of reflex vertex.
!
!    Input, real ( kind = 8 ) ASPC2D, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of
!    a separator.
!
!    Input, real ( kind = 8 ) ATOL2D, the angle tolerance parameter in
!    radians used in accepting separators.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 3 times number of vertices in polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 5 times number of vertices in polygon.
!
!    Input, real ( kind = 8 ) X(1:*), Y(1:*), IANG(1:*), integer LINK(1:*),
!    the data structure for simple polygonal region containing reflex vertex;
!    arrays are for x- and y-coordinates, interior angle, counterclockwise
!    link.
!
!    Output, integer ( kind = 4 ) W1, the index in X,Y of vertex which is the endpoint
!    of separator in inner cone or right cone with respect to reflex vertex;
!    on edge starting at (X(-W1),Y(-W1)) if negative.
!
!    Output, integer ( kind = 4 ) W2, 0 if there is only one separator; else index in X,
!    Y of vertex which is endpoint of 2nd separator in left cone;
!    on edge starting at (X(-W2),Y(-W2)) if negative.
!
!    Output, real ( kind = 8 ) PT1(1:2), coordinates of 1st separator
!    endpoint if W1 < 0.
!
!    Output, real ( kind = 8 ) PT2(1:2), coordinates of 2nd separator
!    endpoint if W2 < 0.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Wokspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angsep
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ivis
  integer ( kind = 4 ) ivor
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) link(*)
  integer ( kind = 4 ) maxn
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvor
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) nvsvrt
  real    ( kind = 8 ) pt1(2)
  real    ( kind = 8 ) pt2(2)
  integer ( kind = 4 ) theta
  integer ( kind = 4 ) v
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) wkang
  real    ( kind = 8 ) x(*)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xvor
  real    ( kind = 8 ) y(*)
  integer ( kind = 4 ) yc
  integer ( kind = 4 ) yvor
!
!  Determine number of vertices in polygon containing reflex vertex.
!
  ierr = 0
  nvrt = 0
  v = vr

  do

    v = link(v)

    if ( v == vr ) then
      exit
    end if

    nvrt = nvrt + 1

  end do

  maxn = nvrt + int(iang(vr)/aspc2d)
!
!  Set up work arrays for routine VISPOL, and determine whether there
!  is enough workspace. XC, YC are d.p. arrays of length NVRT in WK,
!  used for the coordinates of the polygon containing the reflex
!  vertex. MAXN positions are reserved for XC, YC since this is the
!  maximum space required by routine VISVRT. IVIS is an integer array
!  of length MAXN in IWK. IVRT is an integer array of length NVRT in
!  IWK used temporarily for storing indices of vertices in X,Y.
!
  if ( maxiw < maxn + nvrt ) then
    ierr = 6
    return
  else if ( maxwk < maxn + maxn ) then
    ierr = 7
    return
  end if

  ivis = 1
  ivrt = ivis + maxn
  xc = 1
  yc = xc + maxn
  v = link(vr)

  do i = 0, nvrt-1
    wk(xc+i) = x(v)
    wk(yc+i) = y(v)
    iwk(ivrt+i) = v
    v = link(v)
  end do

  call vispol(x(vr),y(vr),nvrt-1,wk(xc),wk(yc),nvis,iwk(ivis), ierr )

  if ( ierr /= 0 ) then
    return
  end if
!
!  XC, YC now contain visibility polygon coordinates. Update MAXN
!  and set up d.p. array THETA of length MAXN in WK for routine
!  VISVRT. Elements of IVIS are changed to indices of X,Y after call.
!
  maxn = maxn - nvrt + nvis + 1
  theta = yc + maxn

  if ( maxwk < theta + maxn - 1 ) then
    ierr = 7
    return
  end if

  call visvrt(aspc2d,x(vr),y(vr),nvis,wk(xc),wk(yc),iwk(ivis), &
    maxn-1,nvsvrt,wk(theta))

  wk(theta+nvsvrt) = iang(vr)

  do i = ivis, ivis+nvsvrt

    v = iwk(i)

    if ( 0 <= v ) then
      iwk(i) = iwk(ivrt+v)
    else
      iwk(i) = -iwk(ivrt-v-1)
    end if

  end do
!
!  XC, YC now contain coord. of visible vertices to be considered
!  as an endpoint of a separator. Set up work arrays for routine
!  VORNBR. Integer array IVOR and d.p. arrays XVOR, YVOR, each of
!  length NVSVRT+1, are added at the end of IWK and WK arrays.
!
  ivor = ivis + nvsvrt + 1
  xvor = theta + nvsvrt + 1
  yvor = xvor + nvsvrt + 1

  if ( maxiw < ivor + nvsvrt ) then
    ierr = 6
    return
  else if ( maxwk < yvor + nvsvrt ) then
    ierr = 7
    return
  end if

  call vornbr(x(vr),y(vr),nvsvrt,wk(xc),wk(yc),nvor,iwk(ivor), &
     wk(xvor),wk(yvor), ierr )

  if ( ierr /= 0 ) then
    return
  end if
!
!  Set up d.p. array WKANG of length NVOR+1 <= NVSVRT+1 in WK for
!  routine FNDSEP. Only Voronoi neighbors are considered as an
!  endpoint of a separator in the first call to FNDSEP. If the
!  minimum angle created at the boundary by the separator(s) is too
!  small, then a second call is made to FNDSEP in which all visible
!  vertices are considered as an endpoint of a separator.
!
  wkang = xvor

  if ( iwk(ivor+nvor) == nvsvrt ) then
    nvor = nvor - 1
  end if

  if ( iwk(ivor) == 0 ) then
    ivor = ivor + 1
    nvor = nvor - 1
  end if

  call fndspf(atol2d+atol2d,x(vr),y(vr),nvsvrt,wk(xc),wk(yc), &
     iwk(ivis),wk(theta),nvor,iwk(ivor),x,y,iang,link,angsep,i1,i2, &
     wk(wkang))

  if ( angsep < atol2d ) then

    ivrt = ivis + nvsvrt + 1

    do i = 1, nvsvrt-1
      iwk(ivrt+i-1) = i
    end do

    call fndspf(atol2d+atol2d,x(vr),y(vr),nvsvrt,wk(xc),wk(yc), &
      iwk(ivis),wk(theta),nvsvrt-2,iwk(ivrt),x,y,iang,link,angsep, &
      i1,i2,wk(wkang))

  end if

  w1 = iwk(ivis+i1)

  if ( w1 < 0 ) then
    pt1(1) = wk(xc+i1)
    pt1(2) = wk(yc+i1)
  end if

  if ( i2 < 0 ) then
    w2 = 0
  else
    w2 = iwk(ivis+i2)
    if ( w2 < 0 ) then
      pt2(1) = wk(xc+i2)
      pt2(2) = wk(yc+i2)
    end if
  end if

  return
end
