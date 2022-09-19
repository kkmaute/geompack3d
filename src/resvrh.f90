subroutine resvrh ( xh, yh, aspc2d, atol2d, nvrt, maxn, maxiw, maxwk, x, y, &
  link, xc, yc, ivrt, v, xv, yv, iwk, wk, ierr )

!*****************************************************************************80
!
!! RESVRH resolves a hole vertex on a face by finding a separator.
!
!  Discussion:
!
!    This routine resolves a top or bottom hole vertex on a face by finding
!    a separator above or below vertex, respectively, where the
!    subregion above or below vertex is given.
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
!    Input, real ( kind = 8 ) XH, YH, the coordinates of hole vertex.
!
!    Input, real ( kind = 8 ) ASPC2D, the angle spacing parameter in
!    radians used in controlling vertices to be considered as an endpoint
!    of a separator.
!
!    Input, real ( kind = 8 ) ATOL2D, the angle tolerance parameter in
!    radians used in accepting separator to resolve a hole on a face.
!
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in subregion containing
!    hole vertex.
!
!    Input, integer ( kind = 4 ) MAXN, at least NVRT+INT(PI/ASPC2D) indicating space
!    available in XC, YC arrays.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about 2*MAXN.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 3*MAXN.
!
!    Input, real ( kind = 8 ) X(1:*), Y(1:*), integer LINK(1:*), used
!    for 2D representation of decomposition of multiply-connected
!    polygonal face of hole polygons; see routine SPDECH.
!
!    Input/output, real ( kind = 8 ) XC(1:MAXN), YC(1:MAXN).  On input,
!    entries 1 through NVRT contain coordinates in counterclockwise order
!    of polygon above or below hole vertex; (XH,YH), (XC(1),YC(1)), and
!    (XC(NVRT),YC(NVRT)) are on horizontal line Y = YH; YC(2), YC(NVRT-1)
!    are both greater than YH (< YH) if top (bottom) vertex.  On output,
!    input values are overwritten; it is assumed that MAXN is sufficiently
!    large.
!
!    Input, integer ( kind = 4 ) IVRT(1:NVRT), the indices in X, Y arrays of XC, YC vertices.
!
!    Output, integer ( kind = 4 ) V, the index in X, Y arrays of second endpoint of
!    separator if positive; on edge starting at (X(-V),Y(-V)) if negative.
!
!    Output, real ( kind = 8 ) XV, YV, the coordinates of second separator
!    endpoint if V < 0.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxn
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) angsep
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivor
  integer ( kind = 4 ) ivrt(nvrt)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ) link(*)
  integer ( kind = 4 ) nmax
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) nvor
  integer ( kind = 4 ) nvsvrt
  integer ( kind = 4 ) v
  real    ( kind = 8 ) wk(maxwk)
  real    ( kind = 8 ) x(*)
  real    ( kind = 8 ) xc(maxn)
  real    ( kind = 8 ) xh
  real    ( kind = 8 ) xv
  integer ( kind = 4 ) xvor
  real    ( kind = 8 ) y(*)
  real    ( kind = 8 ) yc(maxn)
  real    ( kind = 8 ) yh
  real    ( kind = 8 ) yv
  integer ( kind = 4 ) yvor

  ierr = 0

  if ( maxiw < maxn ) then
    ierr = 6
    return
  end if

  call vispol(xh,yh,nvrt-1,xc,yc,nvis,iwk, ierr )

  if ( ierr /= 0 ) then
    return
  end if

  nmax = maxn - nvrt + nvis + 1

  if ( maxwk < nmax ) then
    ierr = 7
    return
  end if
!
!  XC, YC arrays are now overwritten by visibility polygon vertices,
!  then by visible vertices to be considered as separator endpoint.
!  Elements of IVIS = IWK(1:NVSVRT+1) are changed to indices of X, Y.
!  For 1 < I, 0 < IWK(I) = J if ( XC(I),YC(I)) = (X(K),Y(K)) where
!  K = LINK(J), and IWK(I) = J < 0 if ( XC(I),YC(I)) lies in interior
!  of edge starting at (X(-J),Y(-J)). WK(1:NVSVRT+1) contain angles.
!
  call visvrt(aspc2d,xh,yh,nvis,xc,yc,iwk,nmax-1,nvsvrt,wk)

  iwk(1) = ivrt(1)

  do i = 2, nvsvrt+1

    l = iwk(i)

    if ( 0 <= l ) then
      iwk(i) = ivrt(l)
    else
      iwk(i) = -ivrt(-l)
    end if

  end do
!
!  Determine visible vertices which are Voronoi nbrs of (XH,YH).
!
  ivor = nvsvrt + 2
  xvor = nvsvrt + 2
  yvor = xvor + nvsvrt + 1

  if ( maxiw < ivor + nvsvrt ) then
    ierr = 6
    return
  else if ( maxwk < yvor + nvsvrt ) then
    ierr = 7
    return
  end if

  call vornbr(xh,yh,nvsvrt,xc,yc,nvor,iwk(ivor),wk(xvor),wk(yvor), ierr )

  if ( ierr /= 0 ) then
    return
  end if
!
!  Try to find acceptable separator from Voronoi neighbors based on
!  max-min angle criterion. If not successful, find separator by
!  considering all visible vertices as possible endpoint.
!
  if ( iwk(ivor+nvor) == nvsvrt ) then
    nvor = nvor - 1
  end if

  if ( iwk(ivor) == 0 ) then
    ivor = ivor + 1
    nvor = nvor - 1
  end if

  call fndsph(xh,yh,nvsvrt,xc,yc,iwk,wk,nvor,iwk(ivor),x,y,link, &
     angsep,v)

  if ( angsep < atol2d ) then

    iv = nvsvrt + 1

    do i = 1, nvsvrt-1
      iwk(iv+i) = i
    end do

    iv = iv + 1
    call fndsph(xh,yh,nvsvrt,xc,yc,iwk,wk,nvsvrt-2,iwk(iv),x,y, &
      link,angsep,v)

  end if

  v = v + 1

  if ( 0 < iwk(v) ) then
    v = link(iwk(v))
  else
    xv = xc(v)
    yv = yc(v)
    v = iwk(v)
  end if

  return
end
