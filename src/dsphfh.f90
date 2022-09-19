subroutine dsphfh ( aspc2d, atol2d, nvc, nface, nhole, npolh, maxvc, maxfv, &
  maxiw, maxwk, nvert, npf, vcl, facep, factyp, nrml, fvl, eang, hfl, pfl, &
  htsiz, ht, iwk, wk, ierr )

!*****************************************************************************80
!
!! DSPHFH initializes the polyhedral decomposition data structure.
!
!  Discussion:
!
!    This routine initializes the polyhedral decomposition data structure
!    where there may be non-intersecting holes on boundary faces of
!    polyhedral region.  It is assumed head vertex of outer polygon
!    of each face is a strictly convex vertex, and all polygons are
!    simple.  Faces with holes are decomposed into simple polygons.
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
!    Input, real ( kind = 8 ) ASPC2D, angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ATOL2D, angle tolerance parameter in radians
!    used in accepting separator to resolve a hole on a face.
!
!    Input/output, integer ( kind = 4 ) NVC, number of vertex coordinates.
!
!    Input/output, integer ( kind = 4 ) NFACE, number of faces (outer polygon boundaries) in
!    polyhedral decomposition.
!
!    Input, integer ( kind = 4 ) NHOLE, number of holes (inner polygons) on all faces.
!
!    Input, integer ( kind = 4 ) NPOLH, number of polyhedra in decomposition.
!
!    Input, integer ( kind = 4 ) MAXVC, maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXFV, maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXIW, maximum size available for IWK array; should be about
!    max ( 4*(max number of edges in a polyhedron of decomposition),
!    6*max ( NV + 8*NFHOL) ) where NV is number of vertices and
!    NFHOL is number of holes on a multiply-connected face.
!
!    Input, integer ( kind = 4 ) MAXWK, maximum size available for WK array; should be about
!    7*max ( NV + 8*NFHOL ).
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1,1:NFACE+NHOLE+1), the head pointer to vertex
!    indices in FVL for each polygon (face or hole); 1 = FACEP(1,1) < ...
!    < FACEP(1,NFACE+NHOLE+1).
!
!    Input, integer ( kind = 4 ) FACTYP(1:NFACE); face types: useful for specifying types of
!    boundary faces; entries must be greater than or equal to 0.
!
!    Input/output, integer ( kind = 4 ) FACTYP(NFACE+1:NFACE+NHOLE); for I from NFACE+1 to
!    NFACE+NHOLE, FACTYP(I) = +F or -F where 1 <= F <= NFACE is index of
!    face containing hole and sign is + (-) if hole polygon is
!    oriented counterclockwise (CW) in polyhedron when viewed from outside.
!
!    Input, integer ( kind = 4 ) FVL(1,1:*), the vertex indices; those for Ith polygon are in
!    FVL(1,J) for J = FACEP(1,I),...,FACEP(1,I+1)-1; all those
!    for outer polygons must appear before those for holes,
!    and holes must appear in nondecreasing |FACTYP(I)| order.
!
!    Input, integer ( kind = 4 ) HFL(1:NPOLH+1), the head pointer to face indices in PFL
!    for each polyhedron; 1 = HFL(1) < HFL(2) < ... < HFL(NPOLH+1).
!
!    Input, integer ( kind = 4 ) PFL(1,1:*), the signed face indices; those for Ith
!    polyhedron are in PFL(1,J) for J = HFL(I),...,HFL(I+1)-1; the face index
!    must be negated if the ordering of vertices for the face
!    in FVL is in CW order when viewed from outside Ith polyhedron;
!    indices for holes should not be included.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime number
!    which is greater than or equal to NVC+2.
!
!    Output, integer ( kind = 4 ) NVERT, the number of positions used in FVL, EANG arrays.
!
!    Output, integer ( kind = 4 ) NPF, the number of positions used in PFL array.
!
!    Output, integer ( kind = 4 ) FACEP(1:3,1:NFACE); FACEP(1,F) is head pointer to face
!    vertices; FACEP(2,F) and FACEP(3,F) are signed indices of
!    2 polyhedra sharing face F; if F is boundary face then
!    FACEP(3,F) = 0; the sign of the polyhedron index indicates
!    whether face is oriented counterclockwise (positive) or
!    clockwise (negative) in FVL when viewed from outside
!    polyhedron; if interior face, 2 signs are different.
!
!    Output, real ( kind = 8 ) NRML(1:3,1:NFACE), the normals at faces;
!    NRML(*,F) is unit outward normal of face F with its vertices oriented
!    counterclockwise when viewed from outside polyhedron |FACEP(2,F)|.
!
!    Output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list; 6 rows are for
!    LOC, FACN, SUCC, PRED, EDGA, EDGC; first 4 fields are same as that
!    used for convex polyhedron data structure (see routine DSCPH).
!    EDGA and EDGC give information about the edge UV where
!    V = FVL(SUCC,U). Let LU = FVL(LOC,U), LV = FVL(LOC,V),
!    and SF = +1 (-1) if face containing UV in polyhedron P is
!    oriented counterclockwise (CW) when viewed from outside P.
!    Let WX be the edge corresponding to UV in the adjacent face
!    of P, where X = FVL(SUCC,W).  If (LV-LU)*SF greater than 0, then
!    FVL(EDGC,U) = W, FVL(EDGA,W) = U, and EANG(U) is angle at
!    UV between the 2 faces inside P; else FVL(EDGA,U) = W,
!    FVL(EDGC,W) = U, and EANG(W) is the edge angle.  In other
!    words, if P is viewed from outside with edge UV directed
!    upwards from vertex with smaller LOC value to other vertex,
!    then there is a counterclockwise or clockwise rotation in P
!    from face containing UV to other face as indicated by EDGA
!    or EDGC, respectively (A for counterclockwise, C for
!    clockwise). If the counterclockwise or clockwise rotation
!    between 2 faces is exterior to the region, then the EDGA or EDGC
!    value is 0 and EANG value is -1.
!
!    Output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges common
!    to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J) and is
!    determined by EDGC field.
!
!    Output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices for
!    each polyhedron, row 2 used for link; NPF exceeds the input size by at
!    most NHOLE; it is assumed there is enough space.
!
!    Workspace, integer HT(0:HTSIZ-1), the hash table used to find matching
!    occurrences of polyhedron edges by calling routine EDGHT.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) npolh

  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) d
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) edgv
  real    ( kind = 8 ) en(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,nface+nhole+1)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(nface+nhole)
  integer ( kind = 4 ) ff
  logical              fflag
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) g
  logical              gflag
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) hfl(npolh+1)
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) irem
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lc
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) link
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) locfv
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) nfacin
  integer ( kind = 4 ) nfhol
  integer ( kind = 4 ) nfph
  integer ( kind = 4 ) nht
  integer ( kind = 4 ) npf
  real    ( kind = 8 ) nrml(3,nface+nhole)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) sf
  integer ( kind = 4 ) size
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wk(maxwk)
  integer ( kind = 4 ) wrem
  integer ( kind = 4 ) y

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pi2 = 2.0D+00 * pi
  maxedg = maxiw / 4
  nfph = nface + nhole
  nvert = facep(1,nfph+1) - 1
  npf = hfl(npolh+1) - 1
  maxpf = npf + nhole
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  do i = 1, nface
    facep(2,i) = 0
    facep(3,i) = 0
    k = facep(1,i)
    l = facep(1,i+1) - 1
    do j = k, l
      fvl(facn,j) = i
      fvl(succ,j) = j + 1
      fvl(pred,j) = j - 1
      fvl(edga,j) = 0
      fvl(edgc,j) = 0
      eang(j) = -1.0D+00
    end do
    fvl(succ,l) = k
    fvl(pred,k) = l
  end do

  do i = 1, npolh

    k = hfl(i)
    l = hfl(i+1) - 1

    do j = k, l
      pfl(2,j) = j + 1
      f = pfl(1,j)
      p = sign ( i, f )
      f = abs(f)
      if ( facep(2,f) == 0 ) then
         facep(2,f) = p
      else
         facep(3,f) = p
      end if
    end do

    pfl(2,l) = k

  end do

  do f = 1, nface

    if ( 0 < facep(2,f) * facep(3,f) ) then
      ierr = 321
      return
    end if

  end do
!
!  Process holes.
!
  do i = nface+1, nface+nhole-1
    if ( abs(factyp(i+1)) < abs(factyp(i)) ) then
      ierr = 340
      return
    end if
  end do

  do i = nface+1, nface+nhole

    sf = factyp(i)
    f = abs(sf)

    if ( facep(3,f) /= 0 ) then
      ierr = 341
      return
    end if

    fflag = ( facep(2,f) * sf < 0 )
    facep(2,i) = facep(2,f)
    k = facep(1,i)
    l = facep(1,i+1) - 1

    do j = k, l

      fvl(facn,j) = f

      if ( fflag ) then
        fvl(succ,j) = j + 1
        fvl(pred,j) = j - 1
      else
        fvl(succ,j) = j - 1
        fvl(pred,j) = j + 1
      end if

      fvl(edga,j) = 0
      fvl(edgc,j) = 0
      eang(j) = -1.0D+00

    end do

    if ( fflag ) then
      fvl(succ,l) = k
      fvl(pred,k) = l
    else
      fvl(pred,l) = k
      fvl(succ,k) = l
    end if

  end do
!
!  Compute normals for each face from orientation in FACEP(2,*), and
!  check that face vertices (including those on holes) are coplanar.
!
  g = nface + 1

  do f = 1, nface

    if ( 0 < facep(2,f) ) then
      ccw = succ
    else
      ccw = pred
    end if

    j = facep(1,f)
    k = fvl(ccw,j)
    l = fvl(7-ccw,j)
    lb = fvl(loc,j)
    lc = fvl(loc,k)
    la = fvl(loc,l)

    ab(1:3) = vcl(1:3,lb) - vcl(1:3,la)
    ac(1:3) = vcl(1:3,lc) - vcl(1:3,la)

    nrml(1,f) = ab(2) * ac(3) - ab(3) * ac(2)
    nrml(2,f) = ab(3) * ac(1) - ab(1) * ac(3)
    nrml(3,f) = ab(1) * ac(2) - ab(2) * ac(1)

    leng = sqrt(nrml(1,f)**2 + nrml(2,f)**2 + nrml(3,f)**2)

    if ( 0.0D+00 < leng ) then
      nrml(1:3,f) = nrml(1:3,f) / leng
    end if

    d = nrml(1,f)*vcl(1,la)+nrml(2,f)*vcl(2,la)+nrml(3,f)*vcl(3,la)

    if ( abs(d) <= 1.0D+00 ) then
      dtol = tol
    else
      dtol = tol * abs ( d )
    end if

    gflag = .false.

110 continue

    if ( gflag ) then

      if ( nfph < g ) then
        cycle
      end if

      if ( abs(factyp(g)) /= f ) then
        cycle
      end if

      j = facep(1,g)
      l = j
      g = g + 1

    else

      gflag = .true.
      j = fvl(ccw,k)

      if ( j == l) then
        go to 110
      end if

    end if

120 continue

    lb = fvl(loc,j)
    dotp = dot_product ( nrml(1:3,f), vcl(1:3,lb) )

    if ( dtol < abs(dotp - d) ) then
      ierr = 342
      if ( 2 <= msglvl ) then
        write ( *,600) f,lb
      end if
    end if

    j = fvl(ccw,j)

    if ( j /= l ) then
      go to 120
    end if

    go to 110

  end do

  if ( ierr /= 0 ) then
    return
  end if
!
!  Determine EDGA, EDGC fields and compute EANG values.  Temporarily
!  use PFL entries to record hole polygons for polyhedra.
!
  k = npf

  do i = nface+1, nface+nhole
    k = k + 1
    p = abs(facep(2,i))
    j = hfl(p)
    pfl(1,k) = i
    pfl(2,k) = pfl(2,j)
    pfl(2,j) = k
  end do

  do p = 1,npolh

    nht = 0
    i = hfl(p)

150 continue

    if ( i <= npf ) then
      sf = pfl(1,i)
      f = abs(sf)
      ff = f
    else
      ff = pfl(1,i)
      sf = facep(2,ff)
      f = abs(factyp(ff))
    end if

    do j = facep(1,ff),facep(1,ff+1)-1

      la = fvl(loc,j)
      lb = fvl(loc,fvl(succ,j))

      call edght ( la, lb, j, nvc, htsiz, maxedg, hdfree, last, ht, &
        iwk, k, ierr )

      if ( ierr /= 0 ) then
        ierr = 6
        return
      end if

      if ( k <= 0 ) then

        nht = nht + 1

      else

        nht = nht - 1
        g = fvl(facn,k)
        dotp = dot_product ( nrml(1:3,f), nrml(1:3,g) )
        if ( 1.0D+00-tol < abs ( dotp ) ) then
          dotp = sign ( 1.0D+00, dotp )
        end if
        fflag = (abs(facep(2,f)) == p)
        gflag = (abs(facep(2,g)) == p)

        if ( fflag .neqv. gflag ) then
          dotp = -dotp
        end if

        ang = pi - acos(dotp)
!
!  Determine whether edge angle is reflex.
!
        ab(1:3) = vcl(1:3,lb) - vcl(1:3,la)

        en(1) = nrml(2,f) * ab(3) - nrml(3,f) * ab(2)
        en(2) = nrml(3,f) * ab(1) - nrml(1,f) * ab(3)
        en(3) = nrml(1,f) * ab(2) - nrml(2,f) * ab(1)

        if ( fflag .neqv. ( 0 < sf ) ) then
          en(1:3) = -en(1:3)
        end if
!
!  AC = (midpoint of A and B) + EN - A
!
        do l = 1,3
          ac(l) = 0.5D+00*(vcl(l,lb) - vcl(l,la)) + en(l)
        end do

        dotp = ac(1)*nrml(1,g)+ac(2)*nrml(2,g)+ac(3)*nrml(3,g)
        if ( .not. gflag) dotp = -dotp

        if ( 0.0D+00 < dotp ) then
          ang = pi2 - ang
        end if


        if ( 0 < (lb - la)*sf ) then
          fvl(edgc,j) = k
          fvl(edga,k) = j
          eang(j) = ang
        else
          fvl(edga,j) = k
          fvl(edgc,k) = j
          eang(k) = ang
        end if

      end if

    end do

    i = pfl(2,i)

    if ( i /= hfl(p)) then
      go to 150
    end if

    if ( nht /= 0 ) then
      ierr = 322
      return
    end if

  end do

  if ( 0 < nhole ) then
    do p = 1,npolh
      j = hfl(p)
      pfl(2,j) = j + 1
    end do
  end if
!
!  Decompose faces containing holes into simple polygons.
!
  nfacin = nface
  g = nface + 1

210 continue

  if ( nfph < g ) then
    return
  end if

  k = g
  f = abs(factyp(g))
  iwk(1) = facep(1,f)
  j = iwk(1)
  size = 0

  do

    size = size + 1
    j = fvl(succ,j)

    if ( j == iwk(1) ) then
      exit
    end if

  end do

230 continue

  size = size + (facep(1,g+1) - facep(1,g))
  g = g + 1

  if ( g <= nfph ) then
    if ( abs(factyp(g)) == f) then
      go to 230
    end if
  end if

  nfhol = g - k
  size = size + 8*nfhol

  if ( maxiw < 1 + nfhol + 3*size ) then
    ierr = 6
    return
  else if ( maxwk < size + size ) then
    ierr = 7
    return
  end if

  do i = 2,nfhol+1
    iwk(i) = facep(1,k)
    k = k + 1
  end do

  locfv = nfhol + 2
  link = locfv + size
  edgv = link + size
  irem = edgv + size
  y = size + 1
  wrem = y + size

  call spdech(aspc2d,atol2d,nfhol,nvc,nface,nvert,npf,maxvc,nfph, &
    maxfv,maxpf,maxiw-irem+1,maxwk-wrem+1,vcl,facep,factyp,nrml, &
    fvl,eang,hfl,pfl,iwk,wk,wk(y),iwk(locfv),iwk(link), &
    iwk(edgv),iwk(irem),wk(wrem), ierr )

  if ( ierr == 0 ) then
    go to 210
  end if

  600 format (1x,'face ',i5,' contains non-coplanar vertex ',i5)

  return
end
