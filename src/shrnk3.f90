subroutine shrnk3 ( sdist, nface, vcl, hvl, nrml, fvl, eang, maxsv, maxiw, &
  maxwk, nsvc, nsface, nsvert, svcl, shvl, sfvl, iwk, wk, ierr )

!*****************************************************************************80
!
!! SHRNK3 shrinks a convex polyhedron.
!
!  Discussion:
!
!    This routine shrinks a convex polyhedron by distance SDIST, i.e. finds
!    the intersection of half-spaces determined by planes (corresponding
!    to faces) at a distance SDIST inside the polyhedron.
!
!    It is assumed that the vertices of each face are oriented counterclockwise
!    when viewed from outside the polyhedron, and no 2 adjacent faces are
!    coplanar and no 2 adjacent edges are collinear.
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
!    Input, real ( kind = 8 ) SDIST, the shrink distance.
!
!    Input, integer ( kind = 4 ) NFACE, the number of faces in convex polyhedron.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) HVL(1:NFACE), the head vertex list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit outward normals
!    of faces.
!
!    Input, integer ( kind = 4 ) FVL(1:5,1:*), the face vertex list; see routine DSCPH;
!    FVL may contain some unused columns (with <= 0 LOC value).
!
!    Input, real ( kind = 8 ) EANG(1:*), the angles at edges common to 2
!    faces; EANG(I) corresponds to FVL(*,I).
!
!    Input, integer ( kind = 4 ) MAXSV, the maximum size available for SVCL and SFVL
!    arrays; should be greater than or equal to twice number of edges
!    in original or shrunken polyhedron.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be greater than or equal to maximum of MAXSV and
!    NFACE + 2*(max number of edges per face) + 1.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should be
!    greater than or equal to 5*(max number of edges per face + 1).
!
!    Output, integer ( kind = 4 ) NSFACE, the number of faces in shrunken polyhedron
!    (<= NFACE); if NSFACE < 4 then shrunken polyhedron is degenerate or
!    empty.
!
!    Output, integer ( kind = 4 ) NSVC, the size of SVCL array (on output).
!
!    Output, integer ( kind = 4 ) NSVERT, the size of SFVL array.
!
!    Output, real ( kind = 8 ) SVCL(1:3,1:NSVC), SHVL(1:NFACE),
!    SFVL(1:5,1:NSVERT), similar to VCL, HVL, FVL but for shrunken
!    polyhedron (NSFACE must be at least 4); faces of original polyhedron
!    having corresponding faces in shrunken polyhedron are numbered 1 to
!    NSFACE in SHVL (in same increasing order as HVL); indices of faces
!    of original polyhedron not having corresponding faces in shrunken
!    polyhedron are put in SHVL(NSFACE+1:NFACE).
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxsv
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ar
  integer ( kind = 4 ) b
  real    ( kind = 8 ) br
  real    ( kind = 8 ) cmax
  real    ( kind = 8 ) cxy
  real    ( kind = 8 ) cyz
  real    ( kind = 8 ) dr
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ), parameter :: edgv = 5
  logical              empty
  logical              equal
  integer ( kind = 4 ) eshr
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,*)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedge
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) li
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  logical              match
  integer ( kind = 4 ) n
  real    ( kind = 8 ) nrml(3,nface)
  integer ( kind = 4 ) nsface
  integer ( kind = 4 ) nshr
  integer ( kind = 4 ) nsvc
  integer ( kind = 4 ) nsvert
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) r21
  real    ( kind = 8 ) r22
  real    ( kind = 8 ) r31
  real    ( kind = 8 ) r32
  real    ( kind = 8 ) rhs
  real    ( kind = 8 ) sdist
  integer ( kind = 4 ) sfvl(5,maxsv)
  integer ( kind = 4 ) shvl(nface)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) svcl(3,maxsv)
  real    ( kind = 8 ) sxy
  real    ( kind = 8 ) syz
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) xc
  integer ( kind = 4 ) xs
  integer ( kind = 4 ) yc
  real    ( kind = 8 ) yr
  integer ( kind = 4 ) ys
  real    ( kind = 8 ) zr
!
!  For each face, rotate normal vector to (0,0,1). Rotation matrix is
!    [ CXY     -SXY     0   ]
!    [ CYZ*SXY CYZ*CXY -SYZ ]
!    [ SYZ*SXY SYZ*CXY  CYZ ]
!  Rotate face vertices, then apply 2D shrink algorithm,
!  followed by intersection with nonadjacent face-planes.
!  K is index for SFVL and SVCL.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  k = 0
  nsface = 0

  if ( maxiw < nface ) then
    ierr = 6
    return
  end if

  iedge = nface + 1
  eshr = 1

  do f = 1, nface

    if ( abs ( nrml(1,f) ) <= tol ) then
      leng = nrml(2,f)
      cxy = 1.0D+00
      sxy = 0.0D+00
    else
      leng = sqrt(nrml(1,f)**2 + nrml(2,f)**2)
      cxy = nrml(2,f)/leng
      sxy = nrml(1,f)/leng
    end if

    cyz = nrml(3,f)
    syz = leng
    r21 = cyz * sxy
    r22 = cyz * cxy
    r31 = nrml(1,f)
    r32 = nrml(2,f)
    iwk(1:nface) = 1
    iwk(f) = 0
    i = hvl(f)
    n = 0

20  continue

    n = n + 1
    i = fvl(succ,i)

    if ( i /= hvl(f) ) then
      go to 20
    end if

    if ( maxiw < nface + n + n + 1 ) then
      ierr = 6
      return
    else if ( maxwk < 5*n + 4 ) then
      ierr = 7
      return
    end if

    ind = iedge + n + 1
    xc = eshr + n
    yc = xc + n + 1
    xs = yc + n + 1
    ys = xs + n + 1
    j = 0

30  continue

    li = fvl(loc,i)
    wk(xc+j) = cxy*vcl(1,li) - sxy*vcl(2,li)
    wk(yc+j) = r21*vcl(1,li) + r22*vcl(2,li) - syz*vcl(3,li)

    if ( pi - tol <= eang(i) ) then
      wk(eshr+j) = 0.0D+00
    else
      wk(eshr+j) = sdist / tan ( eang(i) * 0.5D+00 )
    end if

    g = fvl(facn,fvl(edgv,i))
    iwk(ind+j) = g
    iwk(g) = 0
    j = j + 1
    i = fvl(succ,i)

    if ( i /= hvl(f) ) then
      go to 30
    end if

    wk(xc+n) = wk(xc)
    wk(yc+n) = wk(yc)
    rhs = nrml(1,f)*vcl(1,li) + nrml(2,f)*vcl(2,li) + &
      nrml(3,f)*vcl(3,li)
    zr = rhs - sdist

    call shrnk2(n,wk(xc),wk(yc),wk(eshr),nshr,wk(xs),wk(ys), &
      iwk(iedge), ierr )

    if ( ierr /= 0 ) then
      return
    end if

    shvl(f) = 0

    if ( nshr == 0 ) then
      cycle
    end if

    if ( maxsv < k + nshr ) then
      ierr = 13
      return
    end if

    nsface = nsface + 1
    i = k + 1
    shvl(f) = i

    do j = 0, nshr-1
      k = k + 1
      svcl(1,k) = wk(xs+j)
      svcl(2,k) = wk(ys+j)
      svcl(3,k) = zr
      sfvl(loc,k) = k
      sfvl(facn,k) = f
      sfvl(succ,k) = k + 1
      sfvl(pred,k) = k - 1
      sfvl(edgv,k) = -iwk(ind+iwk(iedge+j))
    end do

    sfvl(succ,k) = i
    sfvl(pred,i) = k
!
!  Intersect shrunken polygon with other half-spaces.
!
    do j = 1, nface

      if ( iwk(j) /= 1 ) then
        cycle
      end if

      ar = cxy*nrml(1,j) - sxy*nrml(2,j)
      br = r21*nrml(1,j) + r22*nrml(2,j) - syz*nrml(3,j)
      dr = r31*nrml(1,j) + r32*nrml(2,j) + cyz*nrml(3,j)
      li = fvl(loc,hvl(j))
      rhs = nrml(1,j)*vcl(1,li) + nrml(2,j)*vcl(2,li) + &
        nrml(3,j)*vcl(3,li)
      dr = rhs - sdist - dr*zr
!
!  Compute intersection of half-plane AR*X + BR*Y <= DR
!  with convex polygon in plane Z = ZR.
!
      call xpghpl(ar,br,dr,j,maxsv,k,shvl(f),svcl,sfvl,empty, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      if ( empty ) then
        nsface = nsface - 1
        go to 70
      end if

    end do

    i = shvl(f)

    do

      j = sfvl(loc,i)
      yr = svcl(2,j)
      svcl(2,j) = r22*yr - sxy*svcl(1,j) + r32*zr
      svcl(1,j) = cxy*svcl(1,j) + r21*yr + r31*zr
      svcl(3,j) = cyz*zr - syz*yr
      i = sfvl(succ,i)

      if ( i == shvl(f) ) then
        exit
      end if

    end do

70  continue

  end do

  nsvert = k
  nsvc = k

  if ( nsface < 4 ) then
    return
  end if
!
!  Determine EDGV values by searching appropriate faces.
!
  match = .true.

  do f = 1, nface

    i = shvl(f)

    if ( i == 0 ) then
      cycle
    end if

80  continue

    a = sfvl(edgv,i)

    if ( 0 < a ) go to 100

    g = -a
    a = shvl(g)

90  continue

    if ( sfvl(edgv,a) /= -f ) then
    a = sfvl(succ,a)

    if ( a /= shvl(g) ) then
      go to 90
    end if

    match = .false.
    go to 100
    end if
    sfvl(edgv,i) = a
    sfvl(edgv,a) = i
    shvl(g) = sfvl(succ,a)

100 continue

    i = sfvl(succ,i)

    if ( i /= shvl(f)) go to 80

  end do

  if ( .not. match ) then
    ierr = 314
    return
  end if
!
!  Give vertices with same coordinates the same LOC value.
!  Unused (deleted) entries of SFVL have a negative LOC value.
!
  if ( maxiw < nsvert ) then
    ierr = 6
    return
  end if

  do i = 1, nsvert
    iwk(i) = 1
    if ( sfvl(loc,i) < 0 ) then
      iwk(i) = 0
    end if
  end do

  do i = 1, nsvert

    if ( iwk(i) == 0 ) then
      cycle
    end if

    iwk(i) = 0
    li = sfvl(loc,i)
    b = i

130 continue

    b = sfvl(succ,sfvl(edgv,b))

    if ( b == i ) then
      cycle
    end if

    lb = sfvl(loc,b)
    iwk(b) = 0

    equal = .true.

    do j = 1, 3
      cmax = max ( abs ( svcl(j,li) ), abs ( svcl(j,lb) ) )
      if ( tol * cmax < abs(svcl(j,li) - svcl(j,lb)) .and. tol < cmax ) then
        equal = .false.
        exit
      end if
    end do

    if ( equal ) then
      sfvl(loc,b) = li
    else
      sfvl(loc,b) = li
      match = .false.
    end if

    go to 130

  end do
!
!  Keep only referenced vertices of SVCL, update LOC field of SFVL.
!
  iwk(1:nsvc) = 0

  do i = 1, nsvert
    li = sfvl(loc,i)
    if ( 0 < li ) then
      iwk(li) = 1
    end if
  end do

  j = 0

  do i = 1, nsvc
    if ( iwk(i) == 1 ) then
      j = j + 1
      iwk(i) = j
      svcl(1,j) = svcl(1,i)
      svcl(2,j) = svcl(2,i)
      svcl(3,j) = svcl(3,i)
    end if
  end do

  nsvc = j
  do  i = 1,nsvert
    li = sfvl(loc,i)
    if ( 0 < li ) then
      sfvl(loc,i) = iwk(li)
    end if
  end do
!
!  Update SHVL and SFVL(FACN,*) to get consecutive face indices.
!
  if ( nface <= nsface ) then
    return
  end if

  i = 0
  j = 0

  do f = 1,nface

    if ( shvl(f) <= 0 ) then

      i = i + 1
      iwk(i) = f

    else

      j = j + 1

      if ( j < f ) then

        a = shvl(f)
        shvl(j) = a

210     continue

        sfvl(facn,a) = j
        a = sfvl(succ,a)
        if ( a /= shvl(j)) go to 210

      end if

    end if

  end do

  do i = 1, nface-nsface
    shvl(nsface+i) = iwk(i)
  end do

  return
end
