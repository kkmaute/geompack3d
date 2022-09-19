subroutine dspgdc ( nvc, vcl, incr, ncur, nvbc, icur, ivrt, maxhv, maxpv, &
  maxho, npolg, nvert, nhole, nhola, regnum, hvl, pvl, iang, holv, htsiz, &
  maxedg, ht, edge, map, ierr )

!*****************************************************************************80
!
!! DSPGDC initializes the polygonal decomposition data structure.
!
!  Discussion:
!
!    This routine initializes the polygonal decomposition data structure
!    given an initial decomposition of a polygonal region which
!    may have holes and/or cut, separator, and hole interfaces.
!    Holes and hole interfaces must be simple polygons.
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
!    Input, integer ( kind = 4 ) NVC, the number of distinct vertex coordinates in region.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:NVC), the vertex coordinates of
!    boundary curves in arbitrary order.
!
!    Input, integer ( kind = 4 ) INCR, a positive integer greater than or equal to NVC,
!    e.g. 10000, added to some elements of IVRT array.
!
!    Input, integer ( kind = 4 ) NCUR, the number of boundary curves (includes outer boundary
!    curves of subregions and boundary curves of holes
!    and hole interfaces).
!
!    Input, integer ( kind = 4 ) NVBC(1:NCUR), the number of vertices per boundary curve.
!
!    Input, integer ( kind = 4 ) ICUR(1:NCUR), indicates type and location of the curves:
!    ICUR(I) = 0 if Ith curve is outer boundary curve,
!    ICUR(I) = K if Ith curve is a hole and is inside
!    the subregion to the left of Kth curve,
!    ICUR(I) = -K if Ith curve is a hole interface and is
!    inside the subregion to the left of Kth curve.
!    K must be the index of an outer or hole interface
!    boundary curve (hole interfaces may be nested).
!    If the Ith curve is inside more than one subregion
!    due to nesting of hole interfaces, then the subregion
!    to the left of Kth curve must be the smallest
!    subregion containing the Ith curve.
!
!    Input, integer ( kind = 4 ) IVRT(1:NV), indices in VCL of vertices of boundary curves;
!    NV = NVBC(1) + ... + NVBC(NCUR); the vertices of each
!    boundary curve must be in counterclockwise order; the first NVBC(1)
!    positions of IVRT are used for the first curve; the
!    next NVBC(2) positions are used for second curve, etc.
!    If the Ith curve is the outer boundary of a subregion
!    determined from cut and separator interfaces, then the
!    elements of IVRT which correspond to this curve are used
!    both for an index in VCL and indicating the type of the
!    edge joining a vertex and its successor as follows.
!    Let J be in range of positions used for the Ith curve
!    and K be the index in VCL of the coordinates of a vertex
!    of the Ith curve. Consider the edge originating from this
!    vertex. IVRT(J) = -K if the edge is part of a cut or
!    separator interface (i.e. there is a subregion to right
!    of edge). IVRT(J) = K if the edge is part of the outer
!    boundary of the region (i.e. the unbounded exterior of
!    the region is to the right of edge). IVRT(J) = K + INCR
!    if the edge is part of the boundary of a hole (i.e.
!    there is a bounded area to the right of edge which is
!    not in the region. If the Ith curve is the boundary of
!    a hole or hole interface, then only IVRT(J) = K is used.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays,
!    should be greater than or equal to NCUR + (number of hole interfaces).
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays;
!    should be greater than or equal to NVERT (see below).
!
!    Input, integer ( kind = 4 ) MAXHO, the maximum size available for HOLV array; should be
!    greater than or equal to NHOLE*2 + NHOLA (see below).
!
!    Input, integer ( kind = 4 ) HTSIZ, the size of hash table HT; should be a prime
!    number which is about NSC/2 where NSC is number of separator and cut
!    interface edges.
!
!    Input, integer ( kind = 4 ) MAXEDG, the maximum size available for EDGE array;
!    should be at least NSC.
!
!    Output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions, set to
!    number of outer subregion boundaries plus number of hole interfaces.
!
!    Output, integer ( kind = 4 ) NVERT, the number of vertices in PVL, set to NV plus
!    number of vertices in holes and hole interfaces (< 2*NV).
!
!    Output, integer ( kind = 4 ) NHOLE, the number of holes and hole interfaces.
!
!    Output, integer ( kind = 4 ) NHOLA, the number of 'attached' holes; these holes are
!    attached to the outer boundary of a subregion through vertices
!    or cut interfaces and have their edges in consecutive
!    order on the boundary (<= NV/4).
!
!    Output, integer ( kind = 4 ) REGNUM(1:NPOLG), region numbers to left of outer and hole
!    interface boundary curves, which are set to the indices
!    of ICUR or NVBC; this array may be useful in some
!    applications for identifying which original region a
!    subpolygon belongs to.
!
!    Output, integer ( kind = 4 ) HVL(1:NPOLG+NHOLE), the head vertex list; the first NPOLG
!    positions contain the head vertex (index in PVL) of an
!    outer or hole interface boundary curve in which the
!    vertices of the curve are in counterclockwise order in PVL; next
!    NHOLE positions contain the head vertex of a hole or
!    hole interface in which vertices are in CW order in PVL.
!
!    Output, integer ( kind = 4 ) PVL(1:4,1:NVERT), IANG(1:NVERT), the polygon vertex list
!    and interior angles; contains the 5 'arrays' LOC, POLG, SUCC
!    EDGV, IANG (the first 4 are integer arrays, the last
!    is a double precision array); the vertices of each
!    polygon (except for holes) are stored in counterclockwise order in a
!    circular linked list. PVL(LOC,V) is the location in VCL
!    of the coordinates of 'vertex' (index) V. IANG(V) is
!    the interior angle at vertex V. PVL(POLG,V) is polygon
!    number (index of HVL) of subregion containing vertex V
!    (this entry is different from the polygon index only
!    for holes). PVL(SUCC,V) is index in PVL of successor
!    vertex of vertex V. PVL(EDGV,V) gives information about
!    the edge joining vertices V and its successor - if the
!    edge is part of 1 polygon then PVL(EDGV,V) = 0; if the
!    edge is common to 2 polygons then PVL(EDGV,V) greater than 0 and
!    is equal to the index in PVL of the successor vertex
!    as represented in the other polygon; i.e. in latter
!    case, PVL(LOC,PVL(EDGV,V)) = PVL(LOC,PVL(SUCC,V)) and
!    PVL(EDGV,PVL(EDGV,V)) = V.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2+NHOLA), indices in PVL of top or bottom
!    vertex of holes; first (next) NHOLE entries are for top (bottom)
!    vertices of holes and hole interfaces, with top (bottom)
!    vertices sorted in decreasing (increasing) lexicographic
!    (y,x) order of coordinate; last NHOLA entries are for attached
!    holes; if bottom vertex of attached hole is a simple
!    vertex of boundary curve containing the hole then entry
!    contains index of bottom vertex otherwise entry contains
!    index of top vertex (which is simple).
!
!    Workspace, integer MAP(1:NCUR), used for mapping input boundary curve
!    numbers to polygon numbers.
!
!    Workspace, integer HT(0:HTSIZ-1), EDGE(1:4,1:MAXEDG), the hash table
!    and edge records used to determine matching occurrences of separator or
!    cut interface edges by calling routine EDGHT.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg
  integer ( kind = 4 ) maxho
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) ncur
  integer ( kind = 4 ) nvc

  real    ( kind = 8 ) angle
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ), parameter :: edgv = 4
  logical              first
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) hvl(maxhv)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) icur(ncur)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) incr
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) j
  integer ( kind = 4 ) j1
  integer ( kind = 4 ) j2
  integer ( kind = 4 ) jend
  integer ( kind = 4 ) jstr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) kmin
  integer ( kind = 4 ) kpoly
  integer ( kind = 4 ) l
  integer ( kind = 4 ) last
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lvp
  integer ( kind = 4 ) lvs
  integer ( kind = 4 ) map(ncur)
  integer ( kind = 4 ) mpoly
  integer ( kind = 4 ) nh2
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) nhole
  integer ( kind = 4 ) nholi
  integer ( kind = 4 ) nht
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(ncur)
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(maxhv)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) vcl(2,nvc)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) y
  real    ( kind = 8 ) ymax
  real    ( kind = 8 ) ymin

  ierr = 0
  nhola = 0
  nhole = 0
  nholi = 0
  nvert = 0

  do i = 1, ncur
    nvert = nvert + nvbc(i)
    if ( 0 < icur(i) ) then
      nhole = nhole + 1
    else if ( icur(i) < 0 ) then
      nholi = nholi + 1
      nvert = nvert + nvbc(i)
    end if
  end do

  npolg = ncur - nhole
  ipoly = 0
  iv = 0
  nv = 0
  hdfree = 0
  last = 0
  nht = 0
  ht(0:htsiz-1) = 0

  if ( maxhv < ncur + nholi ) then
    ierr = 4
    return
  else if ( maxpv < nvert ) then
    ierr = 5
    return
  else if ( maxho < ( nhole + nholi ) * 2 ) then
    ierr = 2
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for outer boundary curves.
!
  do i = 1, ncur

    if ( icur(i) /= 0 ) then
      map(i) = 0
      go to 40
    end if

    ipoly = ipoly + 1
    regnum(ipoly) = i
    hvl(ipoly) = iv + 1
    map(i) = ipoly
    jstr = nv + 1
    jend = nv + nvbc(i)

    do j = jstr, jend

      iv = iv + 1
      pvl(loc,iv) = abs(ivrt(j))
      pvl(polg,iv) = ipoly
      pvl(succ,iv) = iv + 1

      if ( 0 < ivrt(j) ) then
        pvl(edgv,iv) = 0
      else
!
!  The edge originating from current vertex is on a cut or
!  separator interface. Search in hash table for edge, and
!  insert or delete edge. Set EDGV value if possible.
!
        lv = abs(ivrt(j))

        if ( incr < lv ) then
          lv = lv - incr
        end if

        if ( j < jend ) then
          lvs = abs(ivrt(j+1))
        else
          lvs = abs(ivrt(jstr))
        end if

        if ( incr < lvs ) then
          lvs = lvs - incr
        end if

        call edght ( lv, lvs, iv, nvc, htsiz, maxedg, hdfree, last, ht, &
          edge, ivs, ierr )

        if ( ierr /= 0 ) then
          return
        end if

        if ( 0 < ivs ) then
          pvl(edgv,iv) = ivs
          pvl(edgv,ivs) = iv
          nht = nht - 1
        else
          nht = nht + 1
        end if

      end if

    end do

    pvl(succ,iv) = hvl(ipoly)

40  continue

    nv = nv + nvbc(i)

  end do

  if ( nht /= 0 ) then
    ierr = 215
    return
  end if
!
!  Initialize REGNUM, HVL, PVL arrays for the hole interfaces.
!
  if ( nholi == 0 ) then
    go to 100
  end if

  do i = 1, ncur
    if ( icur(i) < 0 ) then
      ipoly = ipoly + 1
      map(i) = ipoly
    end if
  end do

  nv = 0

  do i = 1, ncur

    if ( 0 < icur(i) ) then
      go to 80
    end if

    ipoly = ipoly + 1
    kpoly = ipoly - nholi
    mpoly = map(-icur(i))
    regnum(kpoly) = i
    hvl(kpoly) = iv + 1
    hvl(ipoly) = iv + 2
    jstr = nv + 1
    jend = nv + nvbc(i)

    do j = jstr, jend
      iv = iv + 2
      pvl(loc,iv-1) = ivrt(j)
      pvl(polg,iv-1) = kpoly
      pvl(succ,iv-1) = iv + 1
      pvl(edgv,iv-1) = iv + 2
      pvl(loc,iv) = ivrt(j)
      pvl(polg,iv) = mpoly
      pvl(succ,iv) = iv - 2
      pvl(edgv,iv) = iv - 3
    end do

    pvl(succ,iv-1) = hvl(kpoly)
    pvl(edgv,iv-1) = hvl(ipoly)
    pvl(succ,hvl(ipoly)) = iv
    pvl(edgv,hvl(ipoly)) = iv - 1

80  continue

    nv = nv + nvbc(i)

  end do
!
!  Initialize HVL, PVL arrays for the ordinary holes.
!
100 continue

  if ( nhole == 0 ) then
    go to 140
  end if

  nv = 0

  do i = 1, ncur

    if ( 0 < icur(i) ) then

      ipoly = ipoly + 1
      mpoly = map(icur(i))
      hvl(ipoly) = iv + 1
      jstr = nv + 1
      jend = nv + nvbc(i)

      do j = jstr, jend
        iv = iv + 1
        pvl(loc,iv) = ivrt(j)
        pvl(polg,iv) = mpoly
        pvl(succ,iv) = iv - 1
        pvl(edgv,iv) = 0
      end do

      pvl(succ,hvl(ipoly)) = iv

    end if

    nv = nv + nvbc(i)

  end do
!
!  Determine bottom or top simple vertex of attached holes.
!
140 continue

  nhole = nhole + nholi
  nh2 = nhole + nhole
  j1 = 0
  j2 = 0

  do i = 1, npolg-nholi

    j = hvl(i)

150 continue

    if ( incr < pvl(loc,j) ) then
      j = pvl(succ,j)
      if ( j /= hvl(i) ) then
        go to 150
      else
        ierr = 216
        return
      end if
    end if

    first = .true.

160 continue

    lv = pvl(loc,j)

    if ( 0 < j1 ) then
      if ( lv <= incr ) then
        j2 = j
      else if ( lv - incr == lvs ) then
        j2 = j
      else
        pvl(loc,j) = lv - incr
      end if
    else if ( incr < lv ) then
      j1 = j
      lvs = lv - incr
      pvl(loc,j) = lvs
    end if

    if ( 0 < j2 ) then
!
!  (Part of) hole starts at vertex J1 and ends at J2.
!
      if ( lv <= incr .and. lv /= lvs ) then
        go to 180
      end if

      k = j1

170   continue

      if ( k == j1 ) then

        kmin = k
        kmax = k
        xmin = vcl(1,lvs)
        ymin = vcl(2,lvs)
        xmax = xmin
        ymax = ymin

      else

        l = pvl(loc,k)
        x = vcl(1,l)
        y = vcl(2,l)

        if ( y < ymin .or. y == ymin .and. x < xmin ) then
          kmin = k
          xmin = x
          ymin = y
        else if ( ymax < y .or. y == ymax .and. xmax < x ) then
          kmax = k
          xmax = x
          ymax = y
        end if

      end if

      k = pvl(succ,k)

      if ( k /= j2) then
        go to 170
      end if

      if ( kmin == j1 ) then
        kmin = kmax
      end if
      nhola = nhola + 1

      if ( maxho < nh2 + nhola ) then
        ierr = 2
        return
      end if

      holv(nh2+nhola) = kmin

180   continue

      j1 = 0
      j2 = 0

      if ( incr < lv ) then
        j1 = j
        pvl(loc,j) = lvs
      end if

    end if

    j = pvl(succ,j)

    if ( first ) then
      first = .false.
      jend = j
      go to 160
    else if ( j /= jend ) then
      go to 160
    end if

  end do
!
!  Initialize IANG array.
!
  do i = 1, npolg+nhole

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

    do

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle(vcl(1,lvp),vcl(2,lvp),vcl(1,lv),vcl(2,lv), &
        vcl(1,lvs),vcl(2,lvs))

      if ( iv == j ) then
        exit
      end if

      lvp = lv
      iv = ivs
      lv = lvs

    end do

  end do
!
!  Initialize HOLV array.
!
  if ( 0 < nhole ) then
    call holvrt(nhole,vcl,hvl(npolg+1),pvl,holv)
  end if

  return
end
