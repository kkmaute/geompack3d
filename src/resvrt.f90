subroutine resvrt ( vr, angspc, angtol, nvc, nvert, maxvc, maxpv, maxiw,  &
  maxwk, vcl, pvl, iang, w1, w2, iwk, wk, ierror )

!*****************************************************************************80
!
!! RESVRT resolves a reflex vertex of a simply connected polygon.
!
!  Discussion:
!
!    This routine resolves a reflex vertex of a simply connected polygon with
!    one or two separators.  The reflex vertex must be a 'simple'
!    vertex of the polygon.
!
!  Modified:
!
!    12 July 1999
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
!    Input, integer ( kind = 4 ) VR, the index in PVL of reflex vertex.
!
!    Input, real ( kind = 8 ) ANGSPC, the angle spacing parameter used in
!    controlling the vertices to be considered as an endpoint of a separator.
!
!    Input, real ( kind = 8 ) ANGTOL, the angle tolerance parameter used in
!    accepting separator(s).
!
!    Input/output, integer ( kind = 4 ) NVC, the number of positions used in VCL array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in PVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 3 times number of vertices in polygon.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 5 times number of vertices in polygon.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) PVL(1:4,1:NVERT), real ( kind = 8 ) IANG(1:NVERT),
!    the polygon vertex list and interior angles.
!
!    Output, integer ( kind = 4 ) W1, the index in PVL of vertex which is the endpoint
!    of separator in inner cone or right cone with respect to the reflex vertex.
!
!    Output, integer ( kind = 4 ) W2, is 0 if there is only one separator; else index
!    in PVL of vertex which is endpoint of second separator in left cone.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) angsep
  real    ( kind = 8 ) angspc
  real    ( kind = 8 ) angtol
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ivis
  integer ( kind = 4 ) ivor
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) maxn
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvor
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) nvsvrt
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) theta
  integer ( kind = 4 ) v
  real    ( kind = 8 ) vcl(2,maxvc)
  integer ( kind = 4 ) vr
  integer ( kind = 4 ) w1
  integer ( kind = 4 ) w2
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) wkang
  integer ( kind = 4 ) xc
  real    ( kind = 8 ) xr
  integer ( kind = 4 ) xvor
  integer ( kind = 4 ) yc
  real    ( kind = 8 ) yr
  integer ( kind = 4 ) yvor
!
!  Determine number of vertices in polygon containing reflex vertex.
!
  nvrt = 0
  v = vr

  do

    v = pvl(succ,v)

    if ( v == vr ) then
      exit
    end if

    nvrt = nvrt + 1

  end do

  maxn = nvrt + int ( iang(vr) / angspc )
  l = pvl(loc,vr)
  xr = vcl(1,l)
  yr = vcl(2,l)
!
!  Set up work arrays for routine VISPOL, and determine whether there
!  is enough workspace. XC, YC are d.p. arrays of length NVRT in WK,
!  used for the coordinates of the polygon containing the reflex
!  vertex. MAXN positions are reserved for XC, YC since this is the
!  maximum space required by routine VISVRT. IVIS is an integer array
!  of length MAXN in IWK. IVRT is an integer array of length NVRT in
!  IWK used temporarily for storing indices of vertices in PVL.
!
  if ( maxiw < maxn + nvrt ) then
    write ( *, * ) ' '
    write ( *, * ) 'RESVRT - Fatal error!'
    write ( *, * ) '  MAXIW < MAXN + NVRT.'
    ierror = 6
    return
  end if

  if ( maxwk < maxn + maxn ) then
    write ( *, * ) ' '
    write ( *, * ) 'RESVRT - Fatal error!'
    write ( *, * ) '  MAXWK < MAXN + MAXN.'
    ierror = 7
    return
  end if

  ivis = 1
  ivrt = ivis + maxn
  xc = 1
  yc = xc + maxn
  v = pvl(succ,vr)

  do i = 0, nvrt-1
    l = pvl(loc,v)
    wk(xc+i) = vcl(1,l)
    wk(yc+i) = vcl(2,l)
    iwk(ivrt+i) = v
    v = pvl(succ,v)
  end do

  call vispol ( xr, yr, nvrt-1, wk(xc), wk(yc), nvis, iwk(ivis), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  XC, YC now contain visibility polygon coordinates. Update MAXN
!  and set up d.p. array THETA of length MAXN in WK for routine
!  VISVRT. Elements of IVIS are changed to indices of PVL after call.
!
  maxn = maxn - nvrt + nvis + 1
  theta = yc + maxn

  if ( maxwk < theta + maxn - 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'RESVRT - Fatal error!'
    write ( *, * ) '  MAXWK < THETA + MAXN - 1.'
    ierror = 7
    return
  end if

  call visvrt ( angspc, xr, yr, nvis, wk(xc), wk(yc), iwk(ivis), maxn-1, &
    nvsvrt, wk(theta) )

  wk(theta+nvsvrt) = iang(vr)

  do i = ivis, ivis+nvsvrt
    l = iwk(i)
    if ( 0 <= l ) then
      iwk(i) = iwk(ivrt+l)
    else
      iwk(i) = -iwk(ivrt-l-1)
    end if
  end do
!
!  XC, YC now contain coordinates of visible vertices to be considered
!  as an endpoint of a separator. Set up work arrays for routine
!  VORNBR. Integer array IVOR and d.p. arrays XVOR, YVOR, each of
!  length NVSVRT+1, are added at the end of IWK and WK arrays.
!
  ivor = ivis + nvsvrt + 1
  xvor = theta + nvsvrt + 1
  yvor = xvor + nvsvrt + 1

  if ( maxiw < ivor + nvsvrt ) then
    ierror = 6
    return
  end if

  if ( maxwk < yvor + nvsvrt ) then
    ierror = 7
    return
  end if

  call vornbr ( xr, yr, nvsvrt, wk(xc), wk(yc), nvor, iwk(ivor), wk(xvor), &
    wk(yvor), ierror )

  if ( ierror /= 0 ) then
    return
  end if
!
!  Set up the array WKANG of length NVOR+1 <= NVSVRT+1 in WK for
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

  call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
    wk(theta), nvor, iwk(ivor), vcl, pvl, iang, angsep, i1, i2, wk(wkang) )

  if ( angsep < angtol ) then

    ivrt = ivis + nvsvrt + 1

    do i = 1, nvsvrt-1
      iwk(ivrt+i-1) = i
    end do

    call fndsep ( angtol+angtol, xr, yr, nvsvrt, wk(xc), wk(yc), iwk(ivis), &
      wk(theta), nvsvrt-2, iwk(ivrt), vcl, pvl, iang, angsep, i1, i2, &
      wk(wkang) )

  end if
!
!  Insert endpoint(s) of separator(s) in vertex coordinate list and
!  polygon vertex list data structures, if they are not yet there.
!
  if ( i2 == -1 ) then
    w2 = 0
  else if ( iwk(ivis+i2) < 0 ) then
    call insvr2 ( wk(xc+i2), wk(yc+i2), -iwk(ivis+i2), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w2, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w2 = iwk(ivis+i2)
  end if

  if ( iwk(ivis+i1) < 0 ) then
    call insvr2 ( wk(xc+i1), wk(yc+i1), -iwk(ivis+i1), nvc, nvert, maxvc, &
      maxpv, vcl, pvl, iang, w1, ierror )
    if ( ierror /= 0 ) then
      return
    end if
  else
    w1 = iwk(ivis+i1)
  end if

  return
end
