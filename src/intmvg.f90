subroutine intmvg ( nsvc, nface, nvert, svcl, hvl, fvl, ibot, itop, h, &
  maxvc, maxwk, nvc, vcl, wk, ierr )

!*****************************************************************************80
!
!! INTMVG generates interior mesh vertices in a shrunken polyhedron.
!
!  Discussion:
!
!    This routine generates interior mesh vertices in a (shrunken) convex
!    polyhedron.  It finds the intersection of the rotated polyhedron (rotated
!    so its diameter is parallel to z-axis) with planes of type z=c
!    at distance H apart.  Then it generates mesh vertices in these
!    convex polygons using a quasi-uniform grid of spacing H.
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
!    Input, integer ( kind = 4 ) NSVC, the number of vertex coordinates for convex
!    polyhedron.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces in convex polyhedron.
!
!    Input, integer ( kind = 4 ) NVERT, the size of FVL array.
!
!    Input/output, real ( kind = 8 ) SVCL(1:3,1:NSVC), the vertex
!    coordinate list.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NFACE), the head vertex list.
!
!    Input, integer ( kind = 4 ) FVL(1:5,1:NVERT), the face vertex list; see routine
!    DSCPH; may contain unused columns, indicated by LOC values <= 0.
!
!    Input, integer ( kind = 4 ) IBOT, ITOP, the indices in SVCL of 2 vertices
!    realizing diameter.
!
!    Input, real ( kind = 8 ) H, the mesh spacing.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be twice maximum number of vertices in intersection of plane with
!    polyhedron.
!
!    Output, integer ( kind = 4 ) NVC, the number of interior mesh vertices generated.
!
!    Output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nsvc
  integer ( kind = 4 ) nvert

  real    ( kind = 8 ) costh
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) dmv(3)
  real    ( kind = 8 ) dz
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,nvert)
  real    ( kind = 8 ) h
  real    ( kind = 8 ) htol
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ibot
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) itop
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) li
  integer ( kind = 4 ) lj
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) n
  integer ( kind = 4 ) np
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcold
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) r
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  real    ( kind = 8 ) r31
  real    ( kind = 8 ) r32
  real    ( kind = 8 ) sinth
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) svcl(3,nsvc)
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) sy
  real    ( kind = 8 ) syz
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wk(maxwk)
  real    ( kind = 8 ) x
  integer ( kind = 4 ) xc
  real    ( kind = 8 ) xl
  real    ( kind = 8 ) xr
  real    ( kind = 8 ) y
  integer ( kind = 4 ) yc
  real    ( kind = 8 ) zr
  real    ( kind = 8 ) zrptol
!
!  Rotate diameter vector to (0,0,1), i.e. rotate vertex coordinates.
!  Rotation matrix is:
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  nvc = 0
  xc = 1
  yc = maxwk / 2 + 1

  dmv(1:3) = svcl(1:3,itop) - svcl(1:3,ibot)

  leng = sqrt ( dmv(1)**2 + dmv(2)**2 + dmv(3)**2 )

  dmv(1:3) = dmv(1:3) / leng

  if ( abs ( dmv(1) ) <= tol ) then
    leng = dmv(2)
    cxy = 1.0D+00
    sxy = 0.0D+00
  else
    leng = sqrt(dmv(1)**2 + dmv(2)**2)
    cxy = dmv(2) / leng
    sxy = dmv(1) / leng
  end if

  cyz = dmv(3)
  syz = leng
  r21 = cyz*sxy
  r22 = cyz*cxy
  r31 = dmv(1)
  r32 = dmv(2)

  do i = 1, nsvc
    x = svcl(1,i)
    y = svcl(2,i)
    svcl(1,i) = cxy*x - sxy*y
    svcl(2,i) = r21*x + r22*y - syz*svcl(3,i)
    svcl(3,i) = r31*x + r32*y + cyz*svcl(3,i)
  end do
!
!  Compute number of planes and intersection of polyhedron with each plane.
!
  np = int ( ( svcl(3,itop) - svcl(3,ibot) ) / h )
  dz = (svcl(3,itop) - svcl(3,ibot) - np*h)*0.5D+00
  htol = h*tol

  if ( dz <= htol ) then
    np = np - 1
    dz = h*0.5D+00
  end if

  if ( np < 0 ) then
    if ( maxvc < 1 ) then
      ierr = 14
    else
      nvc = 1
      x = svcl(1,itop)
      y = svcl(2,itop)
      vcl(1,nvc) = cxy*x + r21*y + r31*svcl(3,itop)
      vcl(2,nvc) = r22*y - sxy*x + r32*svcl(3,itop)
      vcl(3,nvc) = cyz*svcl(3,itop) - syz*y
    end if
    return
  end if

  zr = svcl(3,itop) - dz
  f = 0

  do i = 1, nvert
    if ( fvl(loc,i) == itop ) then
      f = fvl(facn,i)
      hvl(f) = i
      exit
    end if
  end do

  if ( f == 0 ) then
    ierr = 318
    return
  end if

  do k = 0, np

    i = hvl(f)
    li = fvl(loc,i)
    zrptol = zr + htol
!
!  Z-coordinate of vertex I is greater than ZR + HTOL.
!
60  continue

    j = fvl(succ,i)
    lj = fvl(loc,j)

    if ( svcl(3,li) <= svcl(3,lj) ) then
      i = fvl(succ,fvl(edgv,i))
      f = fvl(facn,i)
      go to 60
    else if ( zrptol < svcl(3,lj) ) then
      i = j
      li = lj
      go to 60
    end if
!
!  Trace out convex polygon on plane Z = ZR inside polyhedron.
!
    hvl(f) = i
    nvrt = 0

70  continue

    t = zr - svcl(3,lj)

    if ( t <= htol ) then

      if ( 0 < nvrt .and. svcl(1,lj) == wk(xc) .and. &
        svcl(2,lj) == wk(yc) ) then
        go to 100
      end if

      wk(xc+nvrt) = svcl(1,lj)
      wk(yc+nvrt) = svcl(2,lj)

80    continue

      l = fvl(loc,fvl(succ,j))

      if ( zrptol < svcl(3,l) ) then
        j = fvl(succ,fvl(edgv,j))
        go to 80
      end if

    else

      t = t / (svcl(3,li) - svcl(3,lj))
      wk(xc+nvrt) = svcl(1,lj) + t * (svcl(1,li) - svcl(1,lj))
      wk(yc+nvrt) = svcl(2,lj) + t * (svcl(2,li) - svcl(2,lj))

    end if

    nvrt = nvrt + 1

    if ( yc - 1 <= nvrt ) then
      ierr = 7
      return
    end if

90  continue

    j = fvl(succ,j)
    lj = fvl(loc,j)

    if ( svcl(3,lj) <= zrptol ) then
      go to 90
    end if

    i = fvl(edgv,fvl(pred,j))
    if ( fvl(facn,i) == f ) go to 100
    li = fvl(loc,i)
    j = fvl(succ,i)
    lj = fvl(loc,j)
    go to 70

100 continue

    wk(xc+nvrt) = wk(xc)
    wk(yc+nvrt) = wk(yc)
    call diam2(nvrt,wk(xc+1),wk(yc+1),i1,i2,t,ierr)

    if ( ierr /= 0 ) then
      return
    end if

    call rotpg ( nvrt, wk(xc), wk(yc), i1, i2, ib, costh, sinth )
    n = int((wk(yc) - wk(yc+ib))/h)
    y = wk(yc) - 0.5D+00*(wk(yc) - wk(yc+ib) - real ( n, kind = 8 ) * h )
    l = 0
    r = nvrt
    nvcold = nvc
!
!  Determine left and right x-coordinates of polygon for
!  scan line with y-coordinate Y, and generate mesh vertices.
!
    do i = 0, n

      do while ( y < wk(yc+l+1) )
        l = l + 1
      end do

      do while ( y < wk(yc+r-1) )
        r = r - 1
      end do

      xl = wk(xc+l) + (wk(xc+l+1) - wk(xc+l))*(y - wk(yc+l))/ &
           (wk(yc+l+1) - wk(yc+l))

      xr = wk(xc+r) + (wk(xc+r-1) - wk(xc+r))*(y - wk(yc+r))/ &
           (wk(yc+r-1) - wk(yc+r))

      m = int ( ( xr - xl ) / h )
      x = xl + 0.5D+00*(xr - xl - real ( m, kind = 8 ) * h )
      cy = costh * y
      sy = sinth * y

      if ( maxvc < nvc + m + 1 ) then
        ierr = 14
        return
      end if

      do j = 0, m
        nvc = nvc + 1
        vcl(1,nvc) = costh*x + sy
        vcl(2,nvc) = cy - sinth*x
        x = x + h
      end do

      y = y - h
      if ( y < wk(yc+ib) ) then
        y = wk(yc+ib)
      end if

    end do

    do i = nvcold+1, nvc
      x = vcl(1,i)
      y = vcl(2,i)
      vcl(1,i) = cxy * x + r21 * y + r31 * zr
      vcl(2,i) = r22 * y - sxy * x + r32 * zr
      vcl(3,i) = cyz * zr - syz * y
    end do

    zr = zr - h

  end do

  return
end
