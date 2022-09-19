subroutine trpolg ( nvrt, xc, yc, h, nbc, bndcyc, ldv, nvc, ntri, maxvc,  &
  maxti, maxiw, maxwk, vcl, til, iwk, wk, ierror )

!*****************************************************************************80
!
!! TRPOLG generates a Delaunay triangular mesh inside a convex polygon.
!
!  Discussion:
!
!    This routine generates a Delaunay trianglar mesh inside a convex
!    polygon.  A quasi-uniform grid of spacing H is used.
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
!    Input, integer ( kind = 4 ) NVRT, the number of vertices on the boundary of
!    convex polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates
!    in counterclockwise order; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)); it is
!    assumed that all interior angles are < PI.
!
!    Input, real ( kind = 8 ) H, the spacing of mesh vertices in polygon.
!
!    Input, integer ( kind = 4 ) NBC, the size of BNDCYC.
!
!    Input/output, integer ( kind = 4 ) BNDCYC(0:NBC), the indices in VCL of mesh
!    vertices of boundary cycle; BNDCYC(0) = BNDCYC(NBC); contains
!    (XC(I),YC(I)).
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of coordinates or positions used
!    in VCL array.
!
!    Input/output, integer ( kind = 4 ) NTRI, the number of triangles or positions used
!    in TIL.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXTI, the maximum size available for TIL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array, should
!    be at least 6*(1 + INT(DIAM/H)) + 4*(NBC + NCW) where DIAM is
!    diameter of polygon, NCW is number of edges on boundary
!    of interior triangulation.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array, should
!    be at least 3*NVRT+2.
!
!    Input/output, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxti
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nbc
  integer ( kind = 4 ) nvrt

  integer ( kind = 4 ) bndcyc(0:nbc)
  real    ( kind = 8 ) costh
  integer ( kind = 4 ) cwalk
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) h
  real    ( kind = 8 ) hs
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind
  logical              inter
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) maxcw
  integer ( kind = 4 ) mbc
  integer ( kind = 4 ) ncw
  integer ( kind = 4 ) nshr
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) sdist
  real    ( kind = 8 ) sinth
  real    ( kind = 8 ) smdist
  integer ( kind = 4 ) sptr
  integer ( kind = 4 ) tedg
  integer ( kind = 4 ) til(3,maxti)
  real    ( kind = 8 ) vcl(ldv,maxvc)
  real    ( kind = 8 ) wk(maxwk)
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xi
  integer ( kind = 4 ) xs
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) yi
  real    ( kind = 8 ) yr
  integer ( kind = 4 ) ys

  if ( maxiw < nvrt + 1 ) then
    ierror = 6
    return
  end if

  if ( maxwk < 3*nvrt + 2 ) then
    ierror = 7
    return
  end if

  xs = 1
  ys = xs + nvrt + 1
  sdist = ys + nvrt + 1
  iedge = 1
  hs = h / sqrt ( 2.0D+00 )
  wk(sdist:sdist+nvrt-1) = hs

  call shrnk2 ( nvrt, xc, yc, wk(sdist), nshr, wk(xs), wk(ys), iwk(iedge), &
    ierror )

  if ( ierror /= 0 ) then
    return
  end if

  inter = ( 0 < nshr )

  if ( inter ) then

    call diam2 ( nshr, wk(xs+1), wk(ys+1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    call rotpg ( nshr, wk(xs), wk(ys), i1, i2, ibot, costh, sinth )

    maxcw = 6 * ( 1 + int ( ( wk(ys) - wk(ys+ibot) ) / h ) )

    if ( maxiw < maxcw + 1 ) then
      ierror = 6
      return
    end if

    cwalk = 1

    call inttri ( nshr, wk(xs), wk(ys), h, ibot, costh, sinth, ldv, nvc, ntri, &
      maxvc, maxti, maxcw, vcl, til, ncw, iwk(cwalk), ierror )

    if ( ierror /= 0 ) then
      return
    end if
!
!  Determine the mesh vertex which should be moved to front of
!  BNDCYC - closest to CWALK(0) and also with y-coordinate greater than
!  that of CWALK(0) when rotated if 0 < NCW.
!
    x0 = vcl(1,iwk(cwalk))
    y0 = vcl(2,iwk(cwalk))

    if ( 0 < ncw ) then
      yr = sinth * x0 + costh * y0
    end if

    smdist = 100000.0D+00 * h**2

    do i = 0, nbc-1

      xi = vcl(1,bndcyc(i))
      yi = vcl(2,bndcyc(i))

      if ( 0 < ncw ) then
        if ( sinth * xi + costh * yi <= yr ) then
          cycle
        end if
      end if

      dist = ( xi - x0 )**2 + ( yi - y0 )**2

      if ( dist < smdist ) then
        smdist = dist
        ind = i
      end if

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    nt = nbc + ncw
    tedg = cwalk + ncw + 1

  else

    call diam2 ( nvrt, xc(1), yc(1), i1, i2, dist, ierror )

    if ( ierror /= 0 ) then
      return
    end if

    ind = 0

    do

      if ( nbc <= ind ) then
        exit
      end if

      if ( xc(i1) == vcl(1,bndcyc(ind)) .and. &
           yc(i1) == vcl(2,bndcyc(ind)) ) then
        exit
      end if

      ind = ind + 1

    end do

    call rotiar ( nbc, bndcyc, ind )
    bndcyc(nbc) = bndcyc(0)
    mbc = 1

    do

      if ( nbc <= mbc ) then
        exit
      end if

      if ( xc(i2) == vcl(1,bndcyc(mbc)) .and. &
           yc(i2) == vcl(2,bndcyc(mbc)) ) then
        exit
      end if

      mbc = mbc + 1

    end do

    ind = nbc

    do i = mbc+1, mbc+(nbc-mbc-1)/2
      ind = ind - 1
      i1 = bndcyc(i)
      bndcyc(i) = bndcyc(ind)
      bndcyc(ind) = i1
    end do

    bndcyc(nbc) = bndcyc(mbc)
    nt = nbc - 2
    tedg = 1
!
!  Left boundary chain contains mesh vertices BNDCYC(0:MBC)
!  and right chain contains BNDCYC(0,MBC+1:NBC); MBC < NBC.
!
  end if

  if ( maxti < ntri + nt ) then
    ierror = 9
    return
  else if ( maxiw < tedg + 4*nt - 1 ) then
    ierror = 6
    return
  end if

  if ( inter ) then
    call tmerge ( inter, nbc, ncw, bndcyc, iwk(cwalk), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  else
    call tmerge ( inter, mbc, nbc-mbc, bndcyc, bndcyc(mbc), ldv, vcl, &
      til(1,ntri+1), iwk(tedg), ierror )
  end if

  if ( ierror /= 0 ) then
    return
  end if

  sptr = tedg + 3 * nt

  call cvdtri ( inter, ldv, nt, vcl, til(1,ntri+1), iwk(tedg), iwk(sptr), &
    ierror )

  ntri = ntri + nt

  return
end
