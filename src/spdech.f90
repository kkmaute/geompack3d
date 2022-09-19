subroutine spdech ( aspc2d, atol2d, nfhol, nvc, nface, nvert, npf, maxvc, &
  maxfp, maxfv, maxpf, maxiw, maxwk, vcl, facep, factyp, nrml, fvl, eang, &
  hfl, pfl, headp, x, y, locfv, link, edgv, iwk, wk, ierr )

!*****************************************************************************80
!
!! SPDECH decomposes a face of a polyhedral region.
!
!  Discussion:
!
!    This routine decomposes the face of a polyhedral region into simple
!    polygons, where the face may contain holes, and the outer and inner
!    boundary polygons of the face are simple.
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
!    Input, real ( kind = 8 ) ASPC2D, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ATOL2D, the angle tolerance parameter in
!    radians used in accepting separator to resolve a hole on a face.
!
!    Input, integer ( kind = 4 ) NFHOL, the number of holes on face, must be at least 1.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used in
!    FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPF, the number of positions used in PFL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXFP, the maximum size available for FACEP, FACTYP,
!    NRML arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 3*(NV + 8*NFHOL) where NV is number of vertices on face.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 5*(NV + 8*NFHOL).
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces; outward normal corresponds to counterclockwise
!    traversal of face from polyhedron with index |FACEP(2,F)|.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Input, integer ( kind = 4 ) HFL(1:*), the head pointer to face indices in PFL for
!    each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices
!    for each polyhedron.
!
!    Input, integer ( kind = 4 ) HEADP(0:NFHOL), the first entry is head pointer (index
!    of FVL) of outer polygon of face, other entries are head pointers
!    of hole polygons.
!
!    Workspace, integer HEADP(0:NFHOL), the input values are overwritten.
!
!    Workspace, real ( kind = 8 ) X(1:*), Y(1:*), integer LOCFV(1:*),
!    LINK(1:*), EDGV(1:*), used for 2D representation of multiply-connected
!    polygon; assumed size of each array is at least NV + 8*NFHOL where NV
!    is the number of vertices in multiply-connected polygon.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nfhol

  integer ( kind = 4 ) a
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edgv(*)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) headp(0:nfhol)
  integer ( kind = 4 ) hfl(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ilft
  integer ( kind = 4 ) imax
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) irgt
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) link(*)
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) locfv(*)
  integer ( kind = 4 ) maxn
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfh
  integer ( kind = 4 ) niw
  integer ( kind = 4 ) npf
  real    ( kind = 8 ) nrml(3,maxfp)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  real    ( kind = 8 ) r31
  real    ( kind = 8 ) r32
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) syz
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) w
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) wrem
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xh
  real    ( kind = 8 ) xint
  real    ( kind = 8 ) xlft
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) xrgt
  real    ( kind = 8 ) y(*)
  integer ( kind = 4 ) yc
  real    ( kind = 8 ) yh
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin
  real    ( kind = 8 ) zr

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  if ( nfhol <= 0 ) then
    return
  end if

  if ( maxiw < nfhol ) then
    ierr = 6
    return
  end if
!
!  Rotate normal vector of face to (0,0,1). Rotation matrix applied
!  to face vertices is
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!
  l = fvl(loc,headp(0))
  f = fvl(facn,headp(0))

  if ( 0 < facep(2,f) ) then
    ccw = succ
  else
    ccw = pred
  end if

  if ( tol <= abs ( nrml(1,f) ) ) then
    leng = nrml(2,f)
    cxy = 1.0D+00
    sxy = 0.0D+00
  else
    leng = sqrt ( nrml(1,f)**2 + nrml(2,f)**2 )
    cxy = nrml(2,f) / leng
    sxy = nrml(1,f) / leng
  end if

  cyz = nrml(3,f)
  syz = leng
  r21 = cyz * sxy
  r22 = cyz * cxy
  r31 = nrml(1,f)
  r32 = nrml(2,f)
  zr = r31 * vcl(1,l) + r32 * vcl(2,l) + cyz * vcl(3,l)
  k = 1

  do i = 0, nfhol

    h = headp(i)
    j = h
    headp(i) = k

    do

      l = fvl(loc,j)
      x(k) = cxy*vcl(1,l) - sxy*vcl(2,l)
      y(k) = r21*vcl(1,l) + r22*vcl(2,l) - syz*vcl(3,l)
      locfv(k) = j
      link(k) = k + 1
      edgv(k) = 0
      k = k + 1
      j = fvl(ccw,j)

      if ( j == h ) then
        exit
      end if

    end do

    link(k-1) = headp(i)

  end do

  nv = k - 1
!
!  Determine top and bottom vertices of holes on face, and sort top
!  vertices in decreasing (y,x) order using linear insertion sort.
!
  do i = 1, nfhol

    if ( i == nfhol ) then
      k = nv
    else
      k = headp(i+1) - 1
    end if

    h = headp(i)
    imin = h
    imax = h
    xmin = x(h)
    ymin = y(h)
    xmax = xmin
    ymax = ymin

    do j = h+1, k

      if ( y(j) < ymin .or. ( y(j) == ymin .and. x(j) < xmin ) ) then
        imin = j
        xmin = x(j)
        ymin = y(j)
      else if ( ymax < y(j) .or. ( y(j) == ymax .and. xmax < x(j) ) ) then
        imax = j
        xmax = x(j)
        ymax = y(j)
      end if

    end do

    headp(i) = imax
    iwk(i) = imin

  end do

  do i = 2, nfhol

    h = headp(i)
    xmax = x(h)
    ymax = y(h)
    imin = iwk(i)
    j = i

50  continue

    l = headp(j-1)

    if ( y(l) < ymax .or. ( ymax == y(l) .and. x(l) < xmax ) ) then
      headp(j) = l
      iwk(j) = iwk(j-1)
      j = j - 1
      if ( 1 < j ) then
        go to 50
      end if
    end if

    headp(j) = h
    iwk(j) = imin

  end do
!
!  For each hole, find cut edge from top vertex of hole to a point
!  on outer boundary above top vertex, and update data structures.
!  For each top vertex, find 'closest' vertices on outer boundary
!  which are to left and right of top vertex and on horizontal line
!  through top vertex.  The two closest vertices must be on edges
!  which intersect horizontal line and are partially above line.
!
  inc = int ( pi / aspc2d )
  xlft = 0.0D+00
  xrgt = 0.0D+00

  do i = 1, nfhol

    ilft = 0
    irgt = 0
    h = headp(i)
    xh = x(h)
    yh = y(h)
    j = headp(0)

70  continue

    k = link(j)

    if ( yh < y(k) .and. y(j) <= yh ) then

      if ( y(j) == yh ) then
        xint = x(j)
      else
        xint = (yh - y(j))*(x(k) - x(j))/(y(k) - y(j)) + x(j)
      end if

      if ( xh < xint ) then
        if ( xint < xrgt .or. irgt == 0 ) then
          irgt = j
          xrgt = xint
        end if
      end if

    else if ( yh < y(j) .and. y(k) <= yh ) then

      if ( y(k) == yh ) then
        xint = x(k)
      else
        xint = (yh - y(j))*(x(k) - x(j))/(y(k) - y(j)) + x(j)
      end if

      if ( xint < xh ) then
        if ( xlft < xint .or. ilft == 0 ) then
          ilft = j
          xlft = xint
        end if
      end if

    end if

    j = k

    if ( j /= headp(0) ) then
      go to 70
    end if

    if ( ilft == 0 .or. irgt == 0 ) then
      ierr = 344
      return
    end if

    nvrt = 2
    j = irgt

    do

      j = link(j)
      nvrt = nvrt + 1

      if ( j == ilft ) then
        exit
      end if

    end do

    maxn = nvrt + inc
    niw = nvrt + nfhol

    if ( maxiw < niw ) then
      ierr = 6
      return
    else if ( maxwk < maxn + maxn ) then
      ierr = 7
      return
    end if

    yc = maxn + 1
    wrem = yc + maxn
    wk(1) = xrgt
    wk(maxn+1) = yh
    iwk(nfhol+1) = irgt
    wk(nvrt) = xlft
    wk(maxn+nvrt) = yh
    iwk(niw) = link(ilft)
    j = irgt

    do k = 2,nvrt-1
      j = link(j)
      wk(k) = x(j)
      wk(maxn+k) = y(j)
      iwk(nfhol+k) = j
    end do

    call resvrh(xh,yh,aspc2d,atol2d,nvrt,maxn,maxiw-niw, &
      maxwk-wrem+1,x,y,link,wk,wk(yc),iwk(nfhol+1),v,x(nv+1), &
      y(nv+1),iwk(niw+1),wk(wrem), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( v < 0 ) then

      v = -v
      nv = nv + 1
      link(nv) = link(v)
      link(v) = nv

      if ( maxvc < nvc ) then
        ierr = 14
        return
      end if

      vcl(1,nvc+1) = cxy*x(nv) + r21*y(nv) + r31*zr
      vcl(2,nvc+1) = r22*y(nv) - sxy*x(nv) + r32*zr
      vcl(3,nvc+1) = cyz*zr - syz*y(nv)
      a = locfv(v)
      if ( ccw == pred) a = fvl(pred,a)

      call insvr3(a,nvc,nvert,maxfv,vcl,fvl,eang,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      locfv(nv) = fvl(succ,a)
      w = edgv(v)

      if ( w == 0 ) then

        edgv(nv) = 0
        v = nv

      else

        nv = nv + 1
        x(nv) = x(nv-1)
        y(nv) = y(nv-1)
        link(nv) = link(w)
        link(w) = nv

        if ( locfv(nv-1) == nvert ) then
          locfv(nv) = nvert - 1
        else
          locfv(nv) = nvert
        end if

        edgv(v) = nv
        edgv(nv-1) = w
        edgv(w) = nv - 1
        edgv(nv) = v
        v = nv - 1

      end if

    else

      do j = 1,i-1
        if ( iwk(j) == v ) then
          iwk(j) = -iwk(j)
          exit
        end if
      end do

    end if

    nv = nv + 2
    x(nv-1) = x(h)
    y(nv-1) = y(h)
    link(nv-1) = link(h)
    link(h) = nv
    x(nv) = x(v)
    y(nv) = y(v)
    link(nv) = link(v)
    link(v) = nv - 1
    edgv(nv-1) = edgv(h)
    edgv(h) = v
    edgv(nv) = edgv(v)
    edgv(v) = h

    if ( 0 < edgv(nv-1) ) then
      edgv(edgv(nv-1)) = nv - 1
    end if

    if ( 0 < edgv(nv) ) then
      edgv(edgv(nv)) = nv
    end if

    call inseh3(locfv(h),locfv(v),nvert,maxfv,facep,fvl,eang,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    if ( ccw == succ ) then
      locfv(nv-1) = locfv(h)
      locfv(nv) = locfv(v)
      locfv(h) = nvert - 1
      locfv(v) = nvert
    else
      locfv(nv-1) = nvert -1
      locfv(nv) = nvert
    end if

  end do
!
!  For each hole, find separator from bottom vertex of hole if
!  bottom vertex is not the endpoint of a cut edge ( 0 < IWK(I) ).
!  First sort holes needing separator into increasing (y,x) order.
!  As for cut edges, 2nd endpt of separator is below bottom vertex.
!
  nfh = 0

  do i = 1, nfhol
    if ( 0 < iwk(i) ) then
      nfh = nfh + 1
      iwk(nfh) = iwk(i)
    end if
  end do

  do i = 2, nfh

    h = iwk(i)
    xmin = x(h)
    ymin = y(h)
    j = i

140 continue

    l = iwk(j-1)

    if ( ymin < y(l) .or. ( ymin == y(l) .and. xmin < x(l) ) ) then
      iwk(j) = l
      j = j - 1
      if ( 1 < j ) then
        go to 140
      end if
    end if

    iwk(j) = h

  end do

  do i = 1, nfh

    ilft = 0
    irgt = 0
    h = iwk(i)
    xh = x(h)
    yh = y(h)
    j = link(h)

160 continue

    k = link(j)

    if ( yh <= y(k) .and. y(j) < yh ) then

      if ( y(k) == yh ) then
        xint = x(k)
      else
        xint = (yh - y(j))*(x(k) - x(j))/(y(k) - y(j)) + x(j)
      end if

      if ( xh < xint ) then
        if ( xint < xrgt .or. irgt == 0 ) then
          irgt = j
          xrgt = xint
        end if
      end if

    else if ( yh <= y(j) .and. y(k) < yh ) then

      if ( y(j) == yh ) then
        xint = x(j)
      else
        xint = (yh - y(j))*(x(k) - x(j))/(y(k) - y(j)) + x(j)
      end if

      if ( xint < xh ) then
        if ( xlft < xint .or. ilft == 0 ) then
          ilft = j
          xlft = xint
        end if
      end if

    end if

    j = k
    if ( j /= h) go to 160

    if ( ilft == 0 .or. irgt == 0 ) then
      ierr = 344
      return
    end if

    nvrt = 2
    j = ilft

    do

      j = link(j)
      nvrt = nvrt + 1
      if ( j == irgt ) then
        exit
      end if

    end do

    maxn = nvrt + inc
    niw = nvrt + nfh

    if ( maxiw < niw ) then
      ierr = 6
      return
    else if ( maxwk < maxn + maxn ) then
      ierr = 7
      return
    end if

    yc = maxn + 1
    wrem = yc + maxn
    wk(1) = xlft
    wk(maxn+1) = yh
    iwk(nfh+1) = ilft
    wk(nvrt) = xrgt
    wk(maxn+nvrt) = yh
    iwk(niw) = link(irgt)
    j = ilft

    do k = 2,nvrt-1
      j = link(j)
      wk(k) = x(j)
      wk(maxn+k) = y(j)
      iwk(nfh+k) = j
    end do

    call resvrh(xh,yh,aspc2d,atol2d,nvrt,maxn,maxiw-niw, &
      maxwk-wrem+1,x,y,link,wk,wk(yc),iwk(nfh+1),v,x(nv+1), &
      y(nv+1),iwk(niw+1),wk(wrem), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    if ( v < 0 ) then

      v = -v
      nv = nv + 1
      link(nv) = link(v)
      link(v) = nv

      if ( maxvc <= nvc ) then
        ierr = 14
        return
      end if

      vcl(1,nvc+1) = cxy*x(nv) + r21*y(nv) + r31*zr
      vcl(2,nvc+1) = r22*y(nv) - sxy*x(nv) + r32*zr
      vcl(3,nvc+1) = cyz*zr - syz*y(nv)
      a = locfv(v)
      if ( ccw == pred ) then
        a = fvl(pred,a)
      end if

      call insvr3(a,nvc,nvert,maxfv,vcl,fvl,eang,ierr)

      if ( ierr /= 0 ) then
        return
      end if

      locfv(nv) = fvl(succ,a)
      w = edgv(v)

      if ( w == 0 ) then

        edgv(nv) = 0
        v = nv

      else

        nv = nv + 1
        x(nv) = x(nv-1)
        y(nv) = y(nv-1)
        link(nv) = link(w)
        link(w) = nv

        if ( locfv(nv-1) == nvert ) then
          locfv(nv) = nvert - 1
        else
          locfv(nv) = nvert
        end if

        edgv(v) = nv
        edgv(nv-1) = w
        edgv(w) = nv - 1
        edgv(nv) = v
        v = nv - 1

      end if

    end if

    nv = nv + 2
    x(nv-1) = x(h)
    y(nv-1) = y(h)
    link(nv-1) = link(h)
    link(h) = nv
    x(nv) = x(v)
    y(nv) = y(v)
    link(nv) = link(v)
    link(v) = nv - 1
    edgv(nv-1) = edgv(h)
    edgv(h) = v
    edgv(nv) = edgv(v)
    edgv(v) = h

    if ( 0 < edgv(nv-1) ) then
      edgv(edgv(nv-1)) = nv - 1
    end if

    if ( 0 < edgv(nv) ) then
      edgv(edgv(nv)) = nv
    end if

    a = locfv(h)
    call insed3(a,locfv(v),nface,nvert,npf,maxfp,maxfv,maxpf,facep, &
      factyp,nrml,fvl,eang,hfl,pfl,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    if ( ccw == succ ) then
      locfv(nv-1) = locfv(h)
      locfv(nv) = locfv(v)
      locfv(h) = nvert - 1
      locfv(v) = nvert
    else
      locfv(nv-1) = nvert -1
      locfv(nv) = nvert
    end if

  end do

  return
end
