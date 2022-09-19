subroutine tribfc ( h, nvc, nface, nvert, maxvc, maxbt, maxiw, maxwk, vcl, &
  facep, nrml, fvl, edst, edno, fcst, btst, btl, iwk, wk, ierr )

!*****************************************************************************80
!
!! TRIBFC generates a Delaunay triangulation on polyhedron boundary faces.
!
!  Discussion:
!
!    This routine generates 2D Delaunay triangulations on convex boundary
!    faces of convex polyhedra in polyhedral decomposition data
!    structure according to mesh spacings in H array.
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
!    Input, real ( kind = 8 ) H(1:NPOLH), the mesh spacings for NPOLH
!    polyhedron in data structure.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or positions
!    used in VCL.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces or positions used in FACEP array.
!
!    Input, integer ( kind = 4 ) NVERT, the number of positions used in FVL array.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXBT, the maximum size available for BTL array; should be
!    greater than or equal to the number of triangles generated on
!    boundary faces.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    at least NBC + 1 + (amount needed in routine TRPOLG), where NBC
!    is largest number of mesh vertices on boundary of a face.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    at least 2*NVRT + 2 + (amount needed in routine TRPOLG), where
!    NVRT is largest number of vertices of a face.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate
!    list.  On output, updated by mesh vertices generated on boundary faces.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal vectors
!    for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Output, integer ( kind = 4 ) EDST(1:NVERT), the start location in VCL for mesh
!    vertices on each edge in FVL if there are any, else 0.
!
!    Output, integer ( kind = 4 ) EDNO(1:NVERT), the number of mesh vertices on interior
!    of each edge in FVL; entry is negated if mesh vertices are listed in
!    backward order (according to SUCC) in VCL.
!
!    Output, integer ( kind = 4 ) FCST(1:NFACE+1), the start location in VCL for
!    mesh vertices in interior of each face; last entry indicates where mesh
!    vertices in interior of polyhedra start.
!
!    Output, integer ( kind = 4 ) BTST(1:NFACE+1), the start location in BTL for triangles
!    on each face of FACEP; last entry is (total number of triangles + 1).
!
!    Output, integer ( kind = 4 ) BTL(1:3,1:NBT), the boundary triangle list for triangles
!    generated on all faces of decomposition; NBT = BTST(NFACE+1) - 1;
!    entries are indices of VCL.
!
!    Workspace, integer IWK(1:3,1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK),
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxbt
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) a
  logical              bflag
  integer ( kind = 4 ) btl(3,maxbt)
  integer ( kind = 4 ) btst(nface+1)
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  real    ( kind = 8 ) dz
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edno(nvert)
  integer ( kind = 4 ) edst(nvert)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,nface)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fcst(nface+1)
  integer ( kind = 4 ) fvl(6,nvert)
  real    ( kind = 8 ) h(*)
  real    ( kind = 8 ) hh
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real    ( kind = 8 ) leng
  real    ( kind = 8 ) leng1
  integer ( kind = 4 ) li
  integer ( kind = 4 ) lj
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nbc
  integer ( kind = 4 ) nmv
  real    ( kind = 8 ) nrm1
  real    ( kind = 8 ) nrm2
  real    ( kind = 8 ) nrm3
  real    ( kind = 8 ) nrml(3,nface)
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) nvcb
  integer ( kind = 4 ) p
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  real    ( kind = 8 ) r31
  real    ( kind = 8 ) r32
  integer ( kind = 4 ) s
  integer ( kind = 4 ) start
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) syz
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) u(3)
  real    ( kind = 8 ) v(3)
  real    ( kind = 8 ) vcl(3,nvc)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  real    ( kind = 8 ) xt
  integer ( kind = 4 ) yc
  real    ( kind = 8 ) yt
  real    ( kind = 8 ) zr
!
!  Generate mesh vertices on interior of edges of decomposition.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  edst(1:nvert) = -1

  do i = 1, nvert

    if ( edst(i) /= -1 ) then
      cycle
    end if

    a = i
!
!  Determine whether edge is on the boundary of decomposition.
!
20  continue

    j = fvl(edga,a)

    if ( j == 0 ) then
      bflag = .true.
      f = fvl(facn,a)
      p = abs(facep(2,f))
      hh = h(p)
      n = 1
      go to 30
    else if ( j == i ) then
      bflag = .false.
      a = i
      f = fvl(facn,a)
      p = abs(facep(2,f))
      q = abs(facep(3,f))
      hh = h(p) * h(q)
      n = 2
      go to 30
    end if

    a = j
    go to 20
!
!  Compute mesh spacing for edge from polyhedra containing edge.
!
30  continue

    j = fvl(edgc,a)

40  continue

    k = fvl(edgc,j)

    if ( k /= 0 .and. k /= a ) then

      f = fvl(facn,j)
      r = abs ( facep(2,f) )
      s = abs ( facep(3,f) )

      if ( ( .not. bflag ) .and. n == 2 ) then

        if ( r == p .or. r == q ) then
          p = s
        else
          p = r
        end if

      else

        if ( r == p ) then
          p = s
        else
          p = r
        end if

      end if

      hh = hh * h(p)
      n = n + 1

    end if

    j = k

    if ( j /= 0 .and. j /= a ) then
      go to 40
    end if
!
!  Generate mesh vertices on edge and set EDST, EDNO entries.
!
    hh = hh**( 1.0D+00 / real ( n, kind = 8 ) )
    li = fvl(loc,i)
    lj = fvl(loc,fvl(succ,i))
    dx = vcl(1,lj) - vcl(1,li)
    dy = vcl(2,lj) - vcl(2,li)
    dz = vcl(3,lj) - vcl(3,li)
    leng = sqrt ( dx**2 + dy**2 + dz**2 )
    n = int ( leng / hh )

    if ( real ( n, kind = 8 ) / real ( 2*n+1, kind = 8 ) < leng/hh - n ) then
      n = n + 1
    end if

    n = max ( 0, n-1 )
    start = 0

    if ( 1 <= n ) then

      start = nvc + 1
      dx = dx / real ( n+1, kind = 8 )
      dy = dy / real ( n+1, kind = 8 )
      dz = dz / real ( n+1, kind = 8 )

      if ( maxvc < nvc + n ) then
        ierr = 14
        return
      end if

      do k = 1, n
        nvc = nvc + 1
        vcl(1,nvc) = vcl(1,li) + real ( k, kind = 8 ) * dx
        vcl(2,nvc) = vcl(2,li) + real ( k, kind = 8 ) * dy
        vcl(3,nvc) = vcl(3,li) + real ( k, kind = 8 ) * dz
      end do

    end if

    j = a

!60  continue

    do

      edst(j) = start
      lj = fvl(loc,j)

      if ( lj == li ) then
        edno(j) = n
      else
        edno(j) = -n
      end if

      j = fvl(edgc,j)

      if ( j == 0 .or. j == a ) then
        exit
      end if

    end do

  end do
!
!  Generate Delaunay triangulation on each face of decomposition.
!
  ntri = 0

  do f = 1, nface

    p = abs ( facep(2,f) )
    q = abs ( facep(3,f) )

    if ( q == 0 ) then
      hh = h(p)
    else
      hh = sqrt ( h(p) * h(q) )
    end if

    if ( 0 < facep(2,f) ) then
      nrm1 = nrml(1,f)
      nrm2 = nrml(2,f)
      nrm3 = nrml(3,f)
    else
      nrm1 = -nrml(1,f)
      nrm2 = -nrml(2,f)
      nrm3 = -nrml(3,f)
    end if

    i = facep(1,f)
    li = fvl(loc,i)
    zr = nrm1*vcl(1,li) + nrm2*vcl(2,li) + nrm3*vcl(3,li)
!
!  Equation of face is NRM1*X + NRM2*Y + NRM3*Z = ZR.
!  Rotate normal vector to (0,0,1). Rotation matrix is:
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!
    if ( abs(nrm1) <= tol ) then
      leng = nrm2
      cxy = 1.0D+00
      sxy = 0.0D+00
    else
      leng = sqrt(nrm1**2 + nrm2**2)
      cxy = nrm2 / leng
      sxy = nrm1 / leng
    end if

    cyz = nrm3
    syz = leng
    r21 = cyz * sxy
    r22 = cyz * cxy
    r31 = nrm1
    r32 = nrm2
!
!  Rotate mesh vertices on boundary of face.
!
    n = 0

!80  continue

    do

      n = n + 1
      i = fvl(succ,i)

      if ( i == facep(1,f) ) then
        exit
      end if

    end do

    if ( maxwk < n + n + 2 ) then
      ierr = 7
      return
    end if

    xc = 1
    yc = n + 2
    nvcb = nvc
    n = 0
    i = facep(1,f)
    li = fvl(loc,i)
    lj = fvl(loc,fvl(pred,i))
    u(1) = vcl(1,lj) - vcl(1,li)
    u(2) = vcl(2,lj) - vcl(2,li)
    u(3) = vcl(3,lj) - vcl(3,li)
    leng = u(1)**2 + u(2)**2 + u(3)**2

90  continue

    nvcb = nvcb + 1
    vcl(1,nvcb) = cxy * vcl(1,li) - sxy*vcl(2,li)
    vcl(2,nvcb) = r21 * vcl(1,li) + r22*vcl(2,li) - syz*vcl(3,li)
    vcl(3,nvcb) = li
    j = fvl(succ,i)
    lj = fvl(loc,j)

    v(1:3) = vcl(1:3,lj) - vcl(1:3,li)

    leng1 = v(1)**2 + v(2)**2 + v(3)**2
    dotp = dot_product ( u(1:3), v(1:3) ) / sqrt(leng*leng1)

    if ( -1.0D+00 + tol < dotp ) then
      wk(xc+n) = vcl(1,nvcb)
      wk(yc+n) = vcl(2,nvcb)
      n = n + 1
    end if

    if ( 0 <= edno(i) ) then
      k = edst(i)
      s = 1
    else
      k = edst(i) - edno(i) - 1
      s = -1
    end if

    do r = 1, abs(edno(i))
      nvcb = nvcb + 1
      vcl(1,nvcb) = cxy*vcl(1,k) - sxy*vcl(2,k)
      vcl(2,nvcb) = r21*vcl(1,k) + r22*vcl(2,k) - syz*vcl(3,k)
      vcl(3,nvcb) = k
      k = k + s
    end do

    i = j

    if ( i /= facep(1,f) ) then
      li = lj
      u(1:3) = -v(1:3)
      leng = leng1
      go to 90
    end if

    wk(xc+n) = wk(xc)
    wk(yc+n) = wk(yc)
    nmv = nvcb
    nbc = nvcb - nvc

    if ( maxiw < nbc + 1 ) then
      ierr = 6
      return
    end if

    do i = 1, nbc
      iwk(i) = nvc + i
    end do

    iwk(nbc+1) = iwk(1)
    s = ntri + 1
    btst(f) = s
    fcst(f) = nvc + 1

    call trpolg ( n, wk(xc), wk(yc), hh, nbc, iwk, 3, nmv, ntri, maxvc, maxbt, &
      maxiw-nbc-1, maxwk-yc-n, vcl, btl, iwk(nbc+2), wk(yc+n+1), ierr )

    if ( ierr == 3 ) then
      ierr = 14
    end if

    if ( ierr == 9 ) then
      ierr = 19
    end if

    if ( ierr /= 0 ) then
      return
    end if
!
!  Fix up indices of BTL and rotate interior mesh vertices.
!
    do i = s, ntri
      do j = 1, 3
        r = btl(j,i)
        if ( r <= nvcb ) then
          btl(j,i) = vcl(3,r)
        else
          btl(j,i) = r - nbc
        end if
      end do
    end do

    nmv = nmv - nvcb

    do i = 1, nmv
      xt = vcl(1,nvcb+i)
      yt = vcl(2,nvcb+i)
      vcl(1,nvc+i) = cxy*xt + r21*yt + r31*zr
      vcl(2,nvc+i) = r22*yt - sxy*xt + r32*zr
      vcl(3,nvc+i) = cyz*zr - syz*yt
    end do

    nvc = nvc + nmv

  end do

  btst(nface+1) = ntri + 1
  fcst(nface+1) = nvc + 1

  return
end
