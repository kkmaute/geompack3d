subroutine dsphih ( aspc2d, atol2d, angacc, rdacc, nvc, nface, nvert, &
  npolh, npf, nvch, nfach, ipolh, ifach, holint, maxvc, maxfp, maxfv, &
  maxhf, maxpf, maxiw, maxwk, vcl, facep, factyp, nrml, fvl, eang, hfl, &
  pfl, htsiz, ht, iwk, wk, ierr )

!*****************************************************************************80
!
!! DSPHIH updates the polyhedral decomposition data structure.
!
!  Discussion:
!
!    This routine updates the polyhedral decomposition data structure by
!    adding an interior polyhedron hole to one of the polyhedra,
!    where the hole polyhedron may have holes through it ( 0 <= genus).
!
!    The polyhedron containing the hole is decomposed into
!    2 simple polyhedra by joining the hole to the outer boundary
!    using a cut face; this approach is not guaranteed to work.
!    The interior hole may be a hole interface, i.e. the subregion
!    inside the hole is part of the polyhedral region.
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
!    Input, real ( kind = 8 ) ASPC2D, the angle spacing parameter in radians
!    used in controlling vertices to be considered as an endpoint of a
!    separator.
!
!    Input, real ( kind = 8 ) ATOL2D, the angle tolerance parameter in radians
!    used in accepting separator to resolve a hole on a face.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle in
!    radians produced by a cut face.
!
!    Input, real ( kind = 8 ) RDACC, the minimum acceptable relative distance
!    between a cut plane and vertices not on plane.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates (excluding
!    hole).
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces in decomposition
!    (excluding hole).
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL array
!    (excluding hole).
!
!    Input/output, integer ( kind = 4 ) NPOLH, the number of polyhedra in decomposition.
!
!    Input/output, integer ( kind = 4 ) NPF, the number of positions used in PFL array
!    (excluding hole).
!
!    Input, integer ( kind = 4 ) NVCH, the number of vertex coordinates in hole polyhedron.
!
!    Input, integer ( kind = 4 ) NFACH, the number of faces in hole polyhedron.
!
!    Input, integer ( kind = 4 ) IPOLH, the index of polyhedron containing hole.
!
!    Input, integer ( kind = 4 ) IFACH, the index of face of hole (1 <= IFACH <= NFACH) for
!    which attempt is to be made to find a cut face to join with
!    outer boundary; normal vector of cut face is same as
!    that of hole face; it is assumed that plane containing
!    hole face does not intersect any other part of hole polyhedron
!    (such a face does not exist for all polyhedra); IERR is
!    set to 346 if this cut face is not simple (excluding hole
!    face) or does not meet angle or subedge length criteria.
!
!    Input, logical HOLINT, TRUE iff hole polyhedron is a hole interface.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input, integer ( kind = 4 ) MAXFP, the maximum size available for FACEP, FACTYP,
!    NRML arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXHF, the maximum size available for HFL array.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should be
!    about max ( 4*(number of edges in hole polyh), 6*(NE+NV)), where NE
!    is the number of edges of polyhedron IPOLH intersecting plane
!    through hole face IFACH and NV is the number of vertices on hole face.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about 7*(NE+NV).
!
!    [The following 8 subarrays are as output by routine DSPHDC or
!    DSPHFH, and do not include the hole polyhedron.]
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list: row 1
!    is head pointer, rows 2 and 3 are signed polyhedron indices.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types: useful for
!    specifying types of boundary faces; 0 <= entries.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal vectors
!    for faces; outward normal corresponds to counterclockwise traversal of
!    face from polyhedron with index |FACEP(2,F)|.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges
!    common to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J),
!    determined by EDGC field.
!
!    Input/output, HFL(1:NPOLH), the head pointer to face indices in PFL
!    for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices
!    for each polyhedron; row 2 used for link.
!
!    [The following 5 subarrays are similar to input for DSPHDC
!    for the single hole polyhedron, treated as though the region
!    consists of the hole. Input vertex and face indices in range
!    1 to NVCH and 1 to NFACH, respectively, will be incremented
!    by this routine to avoid conflict with those of polyhedral
!    decomposition. Orientation of faces is also changed.]
!
!    Input, real ( kind = 8 ) VCL(1:3,NVC+1:NVC+NVCH), the vertex coordinate
!    list for hole.
!
!    Input, integer ( kind = 4 ) FACEP(1,NFACE+1:NFACE+NFACH+1), the head pointer to vertex
!    indices in FVL for each hole face; 1 = FACEP(1,NFACE+1)
!    < ... < FACEP(1,NFACE+NFACH+1); head vertex of each face
!    must be a strictly convex vertex.
!
!    Input, FACTYP(NFACE+1:NFACE+NFACH), the face types for faces of hole.
!
!    Input, integer ( kind = 4 ) FVL(1,NVERT+1:*), the vertex indices; those for Ith face
!    of hole are in FVL(1,NVERT+J) for J = FACEP(1,NFACE+I),...,
!    FACEP(1,NFACE+I+1)-1.
!
!    Input, integer ( kind = 4 ) PFL(1,NPF+1:NPF+NFACH), the signed face indices for hole
!    polyhedron; face index must be negated if ordering of vertices for
!    face in FVL is in clockwise order when viewed from outside hole.
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is greater than or equal to NVCH+2.
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
  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxhf
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  real    ( kind = 8 ) ab(3)
  real    ( kind = 8 ) ac(3)
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) aspc2d
  real    ( kind = 8 ) atol2d
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) dir(3,3)
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) en(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  logical              fflag
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) g
  logical              gflag
  integer ( kind = 4 ) headp(0:1)
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) hfhol
  integer ( kind = 4 ) hfint
  integer ( kind = 4 ) hfl(maxhf)
  logical              holint
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iface
  integer ( kind = 4 ) ifach
  integer ( kind = 4 ) ifhol
  integer ( kind = 4 ) ipolh
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) last
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lc
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nfach
  integer ( kind = 4 ) nht
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nv2
  integer ( kind = 4 ) nv3
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvch
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) pt(3)
  real    ( kind = 8 ) rdacc
  logical              rflag
  integer ( kind = 4 ) sf
  integer ( kind = 4 ), parameter :: succ = 3
  integer ( kind = 4 ) tfhol
  integer ( kind = 4 ) tfint
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wk(maxwk)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  pi2 = 2.0D+00 * pi
  maxedg = maxiw / 4
  iface = ifach + nface
  hfhol = npf + 1
  ifhol = npf + ifach
  hdfree = 0
  last = 0
  nht = 0

  ht(0:htsiz-1) = 0

  do i = nface+1, nface+nfach+1
    facep(1,i) = facep(1,i) + nvert
  end do

  do i = nface+1, nface+nfach
    facep(2,i) = 0
    facep(3,i) = 0
    k = facep(1,i)
    l = facep(1,i+1) - 1
    do j = k, l
      fvl(loc,j) = fvl(loc,j) + nvc
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

  nvc = nvc + nvch
  npf = npf + nfach

  do j = hfhol, npf

    if ( 0 < pfl(1,j) ) then
      f = -(pfl(1,j) + nface)
    else
      f = -pfl(1,j) + nface
    end if

    pfl(1,j) = f
    pfl(2,j) = j + 1
    p = sign ( ipolh, f )
    facep(2,abs(f)) = p

  end do

  if ( holint ) then

    if ( maxpf < npf + nfach ) then
      ierr = 17
      return
    end if

    hfint = npf + 1
    tfint = npf + nfach

    do j = hfint, tfint
      pfl(1,j) = -pfl(1,j-nfach)
      pfl(2,j) = j + 1
    end do

    pfl(2,tfint) = hfint

  end if

  nface = nface + nfach
  nvert = facep(1,nface+1) - 1
  tfhol = npf

  if ( ifhol == hfhol ) then
    hfhol = hfhol + 1
  else if ( ifhol == tfhol ) then
    tfhol = tfhol - 1
  else
    pfl(2,ifhol-1) = ifhol + 1
  end if
!
!  Compute normals for each hole face from orientation in FACEP(2,*).
!
  do f = nface-nfach+1, nface

    if ( 0 < facep(2,f) ) then
      ccw = succ
    else
      ccw = pred
    end if

    j = facep(1,f)
    lb = fvl(loc,j)
    lc = fvl(loc,fvl(ccw,j))
    la = fvl(loc,fvl(7-ccw,j))
    ab(1:3) = vcl(1:3,lb) - vcl(1:3,la)
    ac(1:3) = vcl(1:3,lc) - vcl(1:3,la)

    nrml(1,f) = ab(2) * ac(3) - ab(3) * ac(2)
    nrml(2,f) = ab(3) * ac(1) - ab(1) * ac(3)
    nrml(3,f) = ab(1) * ac(2) - ab(2) * ac(1)

    leng = sqrt ( sum ( nrml(1:3,f)**2 ) )

    if ( 0.0D+00 < leng ) then
      nrml(1:3,f) = nrml(1:3,f) / leng
    end if

  end do
!
!  Determine EDGA, EDGC fields + compute EANG values for hole edges.
!
  nht = 0

  do i = npf-nfach+1, npf

    sf = pfl(1,i)
    f = abs(sf)

    do j = facep(1,f),facep(1,f+1)-1

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
        dotp = nrml(1,f)*nrml(1,g) + nrml(2,f)*nrml(2,g) + &
        nrml(3,f)*nrml(3,g)
        if ( 1.0D+00 - tol < abs ( dotp ) ) then
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
        ac(1:3) = 0.5D+00*(vcl(1:3,lb) - vcl(1:3,la)) + en(1:3)

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

  end do

  if ( nht /= 0 ) then
    ierr = 322
    return
  end if
!
!  Determine extreme point of hole face IFACE, and 3 directions on
!  cut plane for routine RESHOL to find starting edge of cut face.
!
  if ( holint ) then
    npf = tfint
  end if

  nrmlc(1:3) = -nrml(1:3,iface)

  j = 1
  if ( abs(nrmlc(1)) < abs(nrmlc(2)) ) then
    j = 2
  end if

  if ( abs(nrmlc(j)) < abs(nrmlc(3)) ) then
    j = 3
  end if

  k = j + 1
  l = j - 1

  if ( 3 < k ) then
    k = 1
  end if

  if ( l < 0 ) then
    l = 3
  end if

  nv = 1
  i = facep(1,iface)
  g = i
  m = fvl(loc,i)
  i = fvl(succ,i)

140 continue

  nv = nv + 1
  j = fvl(loc,i)

  if ( vcl(k,m) < vcl(k,j) .or. &
     ( vcl(k,j) == vcl(k,m) .and. vcl(l,m) < vcl(l,j) ) ) then
    g = i
    m = j
  end if

  i = fvl(succ,i)

  if ( i /= facep(1,iface)) then
    go to 140
  end if

  pt(1:3) = vcl(1:3,m)
  nrmlc(4) = nrmlc(1)*pt(1) + nrmlc(2)*pt(2) + nrmlc(3)*pt(3)
  la = fvl(loc,fvl(succ,g))
  lb = fvl(loc,fvl(pred,g))

  dir(1:3,1) = pt(1:3) - vcl(1:3,la)
  dir(1:3,2) = pt(1:3) - vcl(1:3,lb)
  dir(1:3,3) = dir(1:3,1) + dir(1:3,2)

  do j = 1, 3

    leng = sqrt(dir(1,j)**2 + dir(2,j)**2 + dir(3,j)**2)
    dir(1:3,j) = dir(1:3,j)/leng

    call reshol(ipolh,nrmlc,pt,dir(1,j),angacc,rdacc,nvc,nface, &
      nvert,npolh,npf,maxvc,maxfp,maxfv,maxhf,maxpf,maxiw,maxwk, &
      vcl,facep,factyp,nrml,fvl,eang,hfl,pfl,iwk,wk,rflag, ierr )

    if ( ierr /= 348 ) then
      exit
    end if

  end do

  if ( ierr /= 0 ) then
    return
  end if

  if ( .not. rflag ) then
    ierr = 346
    return
  else if ( maxfv < nvert + nv ) then
    ierr = 15
    return
  end if
!
!  Update PFL entries for hole polyhedron + set data structure for SPDECH.
!
  p = facep(2,iface)
  facep(2,iface) = sign ( npolh, p )
  k = hfl(npolh)
  pfl(2,ifhol) = pfl(2,k)
  pfl(2,k) = ifhol
  k = hfl(ipolh)
  pfl(2,tfhol) = pfl(2,k)
  pfl(2,k) = hfhol

  if ( 0 < p ) then
    ccw = succ
  else
    ccw = pred
  end if

  headp(0) = facep(1,nface)
  headp(1) = nvert + 1
  nvert = nvert + nv
  i = facep(1,iface)

  do j = headp(1), nvert
    fvl(loc,j) = fvl(loc,i)
    fvl(facn,j) = nface
    fvl(succ,j) = j + 1
    fvl(pred,j) = j - 1
    i = fvl(ccw,i)
  end do

  fvl(succ,nvert) = headp(1)
  fvl(pred,headp(1)) = nvert

  if ( ccw == pred ) then
    i = fvl(pred,i)
  end if

  do j = headp(1), nvert

    k = fvl(edgc,i)

    if ( 0 < k ) then
      fvl(edgc,i) = j
      fvl(edga,j) = i
      fvl(edgc,j) = k
      fvl(edga,k) = j
      eang(j) = eang(i) - pi
      eang(i) = pi
    else
      k = fvl(edga,i)
      fvl(edga,i) = j
      fvl(edgc,j) = i
      fvl(edga,j) = k
      fvl(edgc,k) = j
      eang(k) = eang(k) - pi
      eang(j) = pi
    end if

    i = fvl(ccw,i)

  end do

  i = headp(0)

  do

    nv = nv + 1
    i = fvl(succ,i)

    if ( i == headp(0) ) then
      exit
    end if

  end do

  nv = nv + 8
  nv2 = nv + nv
  nv3 = nv2 + nv

  if ( maxiw < nv3 ) then
    ierr = 6
    return
  else if ( maxwk < nv2 ) then
    ierr = 7
    return
  end if

  call spdech(aspc2d,atol2d,1,nvc,nface,nvert,npf,maxvc,maxfp,maxfv, &
    maxpf,maxiw-nv3,maxwk-nv2,vcl,facep,factyp,nrml,fvl,eang,hfl, &
    pfl,headp,wk,wk(nv+1),iwk,iwk(nv+1),iwk(nv2+1),iwk(nv3+1), &
    wk(nv2+1), ierr )

  if ( .not. holint .or. ierr /= 0 ) then
    return
  end if
!
!  Add hole interface to data structure
!
  if ( maxhf <= npolh ) then
    ierr = 18
    return
  end if

  npolh = npolh + 1
  hfl(npolh) = hfint

  do i = hfint, tfint

    f = pfl(1,i)
    facep(3,abs(f)) = sign ( npolh, f )
    f = abs(f)
    j = facep(1,f)

    do

      if ( fvl(edgc,j) == 0 ) then

        k = fvl(edga,j)
        l = fvl(edga,k)

        if ( l == 0 ) then
          fvl(edgc,j) = k
          fvl(edga,k) = j
          eang(j) = pi2 - eang(k)
        else
          fvl(edgc,j) = l
          fvl(edga,l) = j
          eang(j) = pi2 - (eang(k) + eang(l))
        end if

      end if

      j = fvl(succ,j)

      if ( j == facep(1,f) ) then
        exit
      end if

    end do

  end do

  return
end
