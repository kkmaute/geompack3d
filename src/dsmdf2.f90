subroutine dsmdf2 ( hflag, nvc, npolg, maxwk, vcl, hvl, pvl, iang, ivrt, &
  xivrt, widsq, edgval, vrtval, area, wk, ierr )

!*****************************************************************************80
!
!! DSMDF2 sets up a mesh distribution function data structure in 2D.
!
!  Discussion:
!
!    This routine sets up a data structure for a heuristic mesh distribution
!    function from data structure for convex polygon decomposition
!    if HFLAG is TRUE, else set up only IVRT and XIVRT.
!
!    It also computes areas of convex polygons.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
!
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical HFLAG, is TRUE if data structure is to be constructed,
!    FALSE if only IVRT, XIVRT, AREA are to be computed.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates in VCL array.
!
!    Input, integer ( kind = 4 ) NPOLG, the number of polygonal subregions in HVL array.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    2 times maximum number of vertices in any polygon.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NPOLG), the head vertex list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), the polygon vertex list.
!
!    Input, real ( kind = 8 ) IANG(1:*), the interior angles.
!
!    Output, integer ( kind = 4 ) IVRT(1:*), the indices of polygon vertices in VCL,
!    ordered by polygon; same size as PVL.
!
!    Output, integer ( kind = 4 ) XIVRT(1:NPOLG+1), the pointer to first vertex of
!    each polygon in IVRT; vertices of polygon K are IVRT(I) for I from
!    XIVRT(K) to XIVRT(K+1)-1.
!
!    Output, real ( kind = 8 ) WIDSQ(1:NPOLG), the square of width of convex
!    polygons.
!
!    Output, real ( kind = 8 ) EDGVAL(1:*), the value associated with each
!    edge of decomposition; same size as PVL.
!
!    Output, real ( kind = 8 ) VRTVAL(1:NVC), the value associated with each
!    vertex of decomposition.
!
!    [Note: Above 5 arrays are for heuristic mdf data structure.]
!
!    Output, real ( kind = 8 ) AREA(1:NPOLG), the area of convex polygons.
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nvc

  real    ( kind = 8 ) area(npolg)
  real    ( kind = 8 ) areapg
  integer ( kind = 4 ), parameter :: edgv = 4
  real    ( kind = 8 ) edgval(*)
  logical              hflag
  integer ( kind = 4 ) hvl(npolg)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) il
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jl
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nvrt
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pimtol
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,*)
  real    ( kind = 8 ) s
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(2,nvc)
  real    ( kind = 8 ) vrtval(nvc)
  real    ( kind = 8 ) widsq(npolg)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xivrt(npolg+1)
  integer ( kind = 4 ) yc
!
!  Compute area and square of width of polygons.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pimtol = pi - tol

  do k = 1, npolg

    nvrt = 0
    i = hvl(k)

    do

      if ( iang(i) < pimtol ) then
        nvrt = nvrt + 1
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    if ( maxwk < nvrt + nvrt ) then
      ierr = 7
      return
    end if

    xc = 0

    do

      if ( iang(i) < pimtol ) then
        j = pvl(loc,i)
        xc = xc + 1
        wk(xc) = vcl(1,j)
        wk(xc+nvrt) = vcl(2,j)
      end if

      i = pvl(succ,i)

      if ( i == hvl(k) ) then
        exit
      end if

    end do

    xc = 1
    yc = xc + nvrt
    area(k) = areapg(nvrt,wk(xc),wk(yc))*0.5D+00

    if ( hflag ) then
      call width2(nvrt,wk(xc),wk(yc),i,j,widsq(k), ierr )
      if ( ierr /= 0 ) then
        return
      end if
    end if

  end do
!
!  Set up IVRT, XIVRT, EDGVAL, VRTVAL arrays.
!
  l = 1

  do k = 1, npolg

    xivrt(k) = l
    i = hvl(k)
    il = pvl(loc,i)

40  continue

    ivrt(l) = il
    j = pvl(succ,i)
    jl = pvl(loc,j)

    if ( hflag ) then
      s = min ( (vcl(1,jl)-vcl(1,il))**2 + (vcl(2,jl)-vcl(2,il))**2, widsq(k) )
      m = pvl(edgv,i)
      if ( 0 < m ) then
        s = min ( s, widsq(pvl(polg,m)) )
      end if
      edgval(l) = s
    end if

    l = l + 1
    i = j
    il = jl

    if ( i /= hvl(k) ) then
      go to 40
    end if

  end do

  xivrt(npolg+1) = l

  if ( .not. hflag ) then
    return
  end if

  vrtval(1:nvc) = 0.0D+00

  do k = 1, npolg

    j = xivrt(k+1) - 1
    l = j

    do i = xivrt(k), l

      il = ivrt(i)

      if ( vrtval(il) == 0.0D+00 ) then
        vrtval(il) = min ( edgval(i), edgval(j) )
      else
        vrtval(il) = min ( vrtval(il), edgval(i), edgval(j) )
      end if

      j = i

    end do

  end do

  return
end
