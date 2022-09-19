subroutine sfcshp ( p, headp, cntr, angacc, mxcos, dtol, nvc, maxvc, vcl, &
  facep, nrml, fvl, eang, pfl, nrmlc, nce, cedge, cdang, aflag, nv, indv, &
  wvc, ierr )

!*****************************************************************************80
!
!! SFCSHP seeks a separator or cut face in a convex polyhedron.
!
!  Discussion:
!
!    This routine attempts to find a separator or cut face in a convex
!    polyhedron P based on the shape of the polyhedron, one which nearly
!    bisects the diameter of the polyhedron and has a normal close to the
!    diameter vector.
!
!    Accept cut a face if it creates no small angles or short edges.
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
!    Input, integer ( kind = 4 ) P, the polyhedron index.
!
!    Input, integer ( kind = 4 ) HEADP, the head pointer to face of PFL for convex
!    polyhedron P.
!
!    Input, real ( kind = 8 ) CNTR(1:3), the weighted centroid of polyhedron.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) MXCOS, the maximum cosine of angle allowed for
!    angles subtended by new subedges with respect to centroid; exception
!    is that endpoints of start edge may be midpoints of existing edges.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine if
!    point is on cut plane.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC).  On input, vertex
!    coordinate list.  On output, some temporary or permanent entries may
!    be added at end of array.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:*), the edge angles.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for each
!    polyhedron; row 2 used for link.
!
!    Input, integer ( kind = 4 ) NV, the number of vertices in polyhedron.
!
!    Output, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut
!    plane plus right hand side constant term of plane equation (if
!    acceptable).
!
!    Output, integer ( kind = 4 ) NCE, the number of edges in cut face (if acceptable); it is
!    assumed there is enough space in the following two arrays.
!
!    Output, integer ( kind = 4 ) CEDGE(1:2,0:NCE), real ( kind = 8 ) CDANG(1:NCE),
!    information describing cut polygon as output by routine SEPFAC.
!
!    Output, logical AFLAG, TRUE iff separator face is found and acceptable.
!
!    Workspace, integer INDV(1:NV).
!
!    Workspace, real ( kind = 8 ) WVC(1:3,1:NV).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) nv

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa(2)
  logical              aflag
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb(2)
  real    ( kind = 8 ) cdang(*)
  real    ( kind = 8 ) ce
  integer ( kind = 4 ) cedge(2,0:*)
  real    ( kind = 8 ) cn
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) cosacc
  real    ( kind = 8 ) cosmax
  real    ( kind = 8 ) d
  real    ( kind = 8 ) da
  real    ( kind = 8 ) db
  real    ( kind = 8 ) diam
  real    ( kind = 8 ) dinrm(4)
  integer ( kind = 4 ) dir
  real    ( kind = 8 ) distol
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dp
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) eang(*)
  real    ( kind = 8 ) edg(3)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) enr(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fr
  real    ( kind = 8 ), parameter, dimension(3) :: fract = (/ &
    0.5D+00, 0.4D+00, 0.6D+00 /)
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) headp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ia
  integer ( kind = 4 ) ib
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indv(nv)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kvc
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) lbb
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) ll(2)
  integer ( kind = 4 ), parameter :: loc = 1
  real    ( kind = 8 ) mndang
  real    ( kind = 8 ) mpt(3,2)
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mxcos
  integer ( kind = 4 ) nce
  logical              neg
  real    ( kind = 8 ) nrml(3,*)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvcp2
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,*)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) s(2)
  integer ( kind = 4 ) sf
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) v(3)
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wvc(3,nv)
!
!  Find diameter and plane which perpendicularly bisects diameter.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  aflag = .false.
  mndang = min ( pi * 0.9D+00, angacc * 3.0D+00 )
  k = 0
  i = headp

10 continue

  f = abs(pfl(1,i))
  j = facep(1,f)

20 continue

  l = fvl(loc,j)
  jj = 1

30 continue

  if ( jj <= k ) then

    if ( l /= fvl(loc,indv(jj)) ) then
      jj = jj + 1
      go to 30
    end if

  else

    k = jj
    wvc(1,k) = vcl(1,l)
    wvc(2,k) = vcl(2,l)
    wvc(3,k) = vcl(3,l)
    indv(k) = j
    if ( nv <= k ) go to 40

  end if

  j = fvl(succ,j)
  if ( j /= facep(1,f)) go to 20

  i = pfl(2,i)
  if ( i /= headp) go to 10

40 continue

  call diam3 ( nv, wvc, ia, ib, diam )
  diam = sqrt(diam)

  dinrm(1:3) = (wvc(1:3,ib) - wvc(1:3,ia))/diam

  dinrm(4) = 0.5D+00*(dinrm(1)*(wvc(1,ia)+wvc(1,ib)) + dinrm(2)* &
    (wvc(2,ia)+wvc(2,ib)) + dinrm(3)*(wvc(3,ia)+wvc(3,ib)))
  distol = 0.2D+00 * diam
!
!  Attempt to use the following existing edges as a start edge:
!  both endpoints of edge are within distance DISTOL of bisector
!  plane and angle between edge and diameter vector is between
!  60 and 120 degrees.
!
  i = headp

50 continue

  sf = pfl(1,i)
  f = abs(sf)
  a = facep(1,f)
  la = fvl(loc,a)

60 continue

  b = fvl(succ,a)
  lb = fvl(loc,b)
  lbb = lb
  if ( (lb - la)*sf < 0) go to 80
  if ( eang(a) < mndang) go to 80

  da = abs(dinrm(1)*vcl(1,la) + dinrm(2)*vcl(2,la) + &
    dinrm(3)*vcl(3,la) - dinrm(4))
  if ( distol < da ) go to 80

  db = abs(dinrm(1)*vcl(1,lb) + dinrm(2)*vcl(2,lb) + &
    dinrm(3)*vcl(3,lb) - dinrm(4))
  if ( distol < db ) go to 80

  edg(1:3) = vcl(1:3,lb) - vcl(1:3,la)

  leng = sqrt(edg(1)**2 + edg(2)**2 + edg(3)**2)
  dotp = abs((edg(1)*dinrm(1) + edg(2)*dinrm(2) + edg(3)* &
    dinrm(3))/leng)
  if ( 0.5D+00 < dotp ) go to 80

  if ( lb < la ) then
    call i4_swap ( la, lb )
    edg(1:3) = -edg(1:3)
  end if

  cedge(1,0) = lb
  cedge(1,1) = la
  cedge(2,1) = a

  edg(1:3) = edg(1:3) / leng

  fr = fvl(facn,fvl(edgc,a))
  enr(1) = edg(2)*nrml(3,fr) - edg(3)*nrml(2,fr)
  enr(2) = edg(3)*nrml(1,fr) - edg(1)*nrml(3,fr)
  enr(3) = edg(1)*nrml(2,fr) - edg(2)*nrml(1,fr)
  neg = (abs(facep(2,fr)) /= p)

  if ( neg ) then
    enr(1:3) = -enr(1:3)
  end if

  do k = 1, 3

    ang = eang(a) * fract(k)

    if ( msglvl == 4 ) then
      write ( *,600) a,lb,la,fr,p,eang(a)*180.0D+00/ pi, &
        ang * 180.0D+00 / pi
    end if

    cdang(1) = ang
    cn = cos(ang)
    if ( neg) cn = -cn
    ce = sin(ang)

    nrmlc(1:3) = cn*nrml(1:3,fr) + ce*enr(1:3)

    nrmlc(4) = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + nrmlc(3)*vcl(3,la)

    call sepfac(p,cntr,nrmlc,angacc,mxcos,dtol,nvc,maxvc,vcl, &
      facep,nrml,fvl,eang,nce,cedge,cdang,aflag, ierr )

    if ( ierr /= 0 .or. aflag ) then
      return
    end if

  end do

80 continue

  a = b
  la = lbb
  if ( a /= facep(1,f) ) go to 60

  i = pfl(2,i)
  if ( i /= headp ) go to 50
!
!  If above strategy cannot find an acceptable cut face, then try
!  to find a start edge in the interior of a face which has vertices
!  at distance greater than or equal to DISTOL above and below bisector
!  plane; start edge should nearly lie on bisector plane.
!
  cosacc = max ( mxcos, cos ( angacc ) )

  if ( maxvc < nvc + 4 ) then
    ierr = 14
    return
  end if

  nvcp2 = nvc + 2
  i = headp

90 continue

  sf = pfl(1,i)
  f = abs(sf)
  j = facep(1,f)
  l = fvl(loc,j)
  da = dinrm(1)*vcl(1,l) + dinrm(2)*vcl(2,l) + dinrm(3)*vcl(3,l) &
    - dinrm(4)
  db = da
  ia = j
  ib = j
  j = fvl(succ,j)

100 continue

  l = fvl(loc,j)
  d = dinrm(1)*vcl(1,l) + dinrm(2)*vcl(2,l) + dinrm(3)*vcl(3,l) - dinrm(4)

  if ( da < d ) then
    da = d
    ia = j
  else if ( d < db ) then
    db = d
    ib = j
  end if

  j = fvl(succ,j)

  if ( j /= facep(1,f) ) then
    go to 100
  end if

  if ( da < distol .or. -distol < db ) go to 160

  dir = succ

  do k = 1, 2

    if ( k == 2 ) then
      dir = pred
    end if

    dp = da
    j = fvl(dir,ia)

110 continue

    l = fvl(loc,j)
    d = dinrm(1)*vcl(1,l) + dinrm(2)*vcl(2,l) + dinrm(3)*vcl(3,l) - dinrm(4)

    if ( 0.0D+00 <= d ) then
      dp = d
      j = fvl(dir,j)
      go to 110
    end if

    bb(k) = j
    aa(k) = fvl(7-dir,j)

    if ( distol <= dp - d ) then

      la = fvl(loc,aa(k))

      mpt(1:3,k) = 0.5D+00 * ( vcl(1:3,la) + vcl(1:3,l) )

      if ( 0.0D+00 <= 0.5D+00*(dp + d) ) then
        aa(k) = -aa(k)
      else
        bb(k) = -bb(k)
      end if

    end if

  end do

  do j = 1,2

    if ( j == 1 ) then

      if ( aa(1) == ia .or. bb(2) == ib ) then
        cycle
      end if

      s(1) = aa(1)
      s(2) = bb(2)

    else

      if ( bb(1) == ib .or. aa(2) == ia ) then
        cycle
      end if

      s(1) = bb(1)
      s(2) = aa(2)

    end if

    if ( 0 < s(1) ) then
      ll(1) = fvl(loc,s(1))
    else
      ll(1) = nvcp2 + 1
      vcl(1:3,ll(1)) = mpt(1:3,1)
    end if

    if ( 0 < s(2) ) then
      ll(2) = fvl(loc,s(2))
    else
      ll(2) = nvcp2 + 2
      vcl(1:3,ll(2)) = mpt(1:3,2)
    end if

    edg(1:3) = vcl(1:3,ll(2)) - vcl(1:3,ll(1))

    leng = edg(1)**2 + edg(2)**2 + edg(3)**2
    cosmax = -1.0D+00

    do k = 1,2

      if ( 0 < s(k) ) then
        l = fvl(loc,fvl(succ,s(k)))
      else
        l = fvl(loc,-s(k))
      end if

      v(1:3) = vcl(1:3,l) - vcl(1:3,ll(k))

      d = v(1)**2 + v(2)**2 + v(3)**2
      dotp = (edg(1)*v(1)+edg(2)*v(2)+edg(3)*v(3))/sqrt(leng*d)
      if ( k == 2) dotp = -dotp
      cosmax = max ( cosmax, dotp )

      if ( s(k) < 0 ) then
        cosmax = max ( cosmax, -dotp )
      else
        l = fvl(loc,fvl(pred,s(k)))

        v(1:3) = vcl(1:3,l) - vcl(1:3,ll(k))

        d = v(1)**2 + v(2)**2 + v(3)**2
        dotp = (edg(1)*v(1) + edg(2)*v(2) + edg(3)*v(3))/ sqrt(leng*d)

        if ( k == 2 ) then
          dotp = -dotp
        end if

        cosmax = max ( cosmax, dotp )

      end if

    end do

    if ( cosacc < cosmax ) then
      cycle
    end if

    neg = ( abs(facep(2,f) ) /= p )
    kvc = nvc

    if ( ll(2) <= nvc ) then

      lb = ll(2)

    else

      kvc = kvc + 1
      lb = kvc

      vcl(1:3,kvc) = vcl(1:3,ll(2))

      b = abs(bb(2))
      k = fvl(loc,b)
      l = fvl(loc,fvl(succ,b))

      if ( 0 < (l - k)*sf ) then
        b = fvl(edgc,b)
      else
        b = fvl(edga,b)
      end if

      cedge(2,0) = b

    end if

    if ( ll(1) <= nvc ) then
      la = ll(1)
      a = s(1)
      if ( 0 < sf ) then
        a = fvl(pred,a)
      end if
    else
      kvc = kvc + 1
      la = kvc

      vcl(1:3,kvc) = vcl(1:3,ll(1))

      a = abs(aa(1))
    end if

    cedge(1,0) = lb
    cedge(1,1) = la
    cedge(2,1) = -a
    leng = sqrt(leng)
    edg(1:3) = edg(1:3)/leng
    enr(1) = edg(2)*nrml(3,f) - edg(3)*nrml(2,f)
    enr(2) = edg(3)*nrml(1,f) - edg(1)*nrml(3,f)
    enr(3) = edg(1)*nrml(2,f) - edg(2)*nrml(1,f)

    if ( neg ) then
      enr(1:3) = -enr(1:3)
    end if

    do k = 1,3

      ang = pi * fract(k)

      if ( msglvl == 4 ) then
        write ( *,610) a,lb,la,f,p, ang*180.0D+00/ pi
      end if

      cdang(1) = ang
      cn = cos(ang)
      if ( neg) cn = -cn
      ce = sin(ang)

      nrmlc(1:3) = cn * nrml(1:3,f) + ce * enr(1:3)

      nrmlc(4) = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + &
        nrmlc(3)*vcl(3,la)

      call sepfac(p,cntr,nrmlc,angacc,mxcos,dtol,nvc,maxvc,vcl, &
        facep,nrml,fvl,eang,nce,cedge,cdang,aflag, ierr )

      if ( ierr /= 0 .or. aflag ) then
        return
      end if

    end do

  end do

160 continue

  i = pfl(2,i)
  if ( i /= headp) go to 90

  600 format (' sfcshp: a,lb,la,f,p,eang(a),ang =',4i5,i6,2f9.4)
  610 format (' sfcshp: a,lb,la,f,p,ang =',5i5,f9.4)

  return
end
