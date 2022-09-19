subroutine reshol ( p, nrmlc, pt, dir, angacc, rdacc, nvc, nface, nvert, &
  npolh, npf, maxvc, maxfp, maxfv, maxhf, maxpf, maxiw, maxwk, vcl, facep, &
  factyp, nrml, fvl, eang, hfl, pfl, iwk, wk, rflag, ierr )

!*****************************************************************************80
!
!! RESHOL finds a cut face in a polyhedron to resolve an interior hole.
!
!  Discussion:
!
!    This routine finds a cut face in polyhedron P to resolve an interior
!    polyhedral hole where a cut plane contains an extreme face of the hole.
!
!    Accept a cut face if it creates no small angles, it is an outer
!    boundary, there are no holes (inner polygons) in it, etc.
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
!    Input, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut
!    plane plus right hand side constant term of plane equation.
!
!    Input, real ( kind = 8 ) PT(1:3), the point (of hole) on cut plane.
!
!    Input, real ( kind = 8 ) DIR(1:3), the unit direction of ray on cut
!    plane from PT for which first intersection with boundary is used to
!    start cut face.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) RDACC, the minimum acceptable relative
!    distance between a cut plane and vertices not on plane.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used
!    in FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions used
!    in HFL array.
!
!    Input/output, integer ( kind = 4 ) NPF, the number of positions used in PFL array.
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
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array; should
!    be about max ( 5*NE, 3*NE+NV ) where NE is number of edges of polyhedron
!    intersecting cut plane and NV is maximum number of face vertices.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about NE.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face indices for
!    each polyhedron.
!
!    Output, logical RFLAG, TRUE iff satisfactory cut polygon is found.
!
!    Workspace, integer IWK(1:MAXIW).
!
!    Workspace, real ( kind = 8 ) WK(1:MAXWK).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxhf
  integer ( kind = 4 ) maxiw
  integer ( kind = 4 ) maxpf
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) maxwk

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angn
  integer ( kind = 4 ) b
  integer ( kind = 4 ) ca
  integer ( kind = 4 ) cb
  integer ( kind = 4 ) ccw
  real    ( kind = 8 ) cmax
  real    ( kind = 8 ) cmin
  real    ( kind = 8 ) cp(3)
  real    ( kind = 8 ) da
  real    ( kind = 8 ) db
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dtol
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fmin
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) hfl(maxhf)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ibeg
  integer ( kind = 4 ) iend
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ineg
  integer ( kind = 4 ) inout
  integer ( kind = 4 ) iomin
  integer ( kind = 4 ) ipos
  real    ( kind = 8 ) ipt(3)
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l(0:1)
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  real    ( kind = 8 ) leng
  integer ( kind = 4 ) link
  integer ( kind = 4 ), parameter :: loc = 1
  real    ( kind = 8 ) mindis
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  integer ( kind = 4 ) ne
  integer ( kind = 4 ) nedg
  integer ( kind = 4 ) nep1
  integer ( kind = 4 ) nface
  logical              nflag
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) pt(3)
  integer ( kind = 4 ) ptr
  real    ( kind = 8 ) ray(3)
  real    ( kind = 8 ) rdacc
  logical              rflag
  integer ( kind = 4 ) sf
  real    ( kind = 8 ) smin
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tmin
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) va(3)
  real    ( kind = 8 ) vb(3)
  integer ( kind = 4 ) vbeg
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) vend
  real    ( kind = 8 ) wk(maxwk)
!
!  Determine edges in polyhedron P and average edge length.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  rflag = .false.
  nedg = 0
  sum2 = 0.0D+00
  ptr = hfl(p)

10 continue

  f = pfl(1,ptr)

  if ( 0 < f ) then
    ccw = succ
  else
    ccw = pred
    f = -f
  end if

  a = facep(1,f)
  lb = fvl(loc,a)

20 continue

  la = lb
  b = fvl(ccw,a)
  lb = fvl(loc,b)

  if ( la < lb ) then
    leng = sqrt((vcl(1,la) - vcl(1,lb))**2 + (vcl(2,la) - &
      vcl(2,lb))**2 + (vcl(3,la) - vcl(3,lb))**2)
    sum2 = sum2 + leng
    nedg = nedg + 1
  end if

  a = b
  if ( a /= facep(1,f)) then
    go to 20
  end if

  ptr = pfl(2,ptr)

  if ( ptr /= hfl(p)) then
    go to 10
  end if

  sum2 = sum2 / real ( nedg, kind = 8 )
  mindis = rdacc * sum2
  dtol = tol * sum2
!
!  Determine which edges of polyhedron intersect cut plane;
!  reject plane if there is a vertex within distance MINDIS of
!  plane which does not lie on plane.
!
  nedg = 0
  ne = 0
  nv = 0
  ptr = hfl(p)

30 continue

  f = pfl(1,ptr)

  if ( 0 < f ) then
    ccw = succ
  else
    ccw = pred
    f = -f
  end if

  a = facep(1,f)
  lb = fvl(loc,a)
  k = 0

40 continue

  la = lb
  b = fvl(ccw,a)
  lb = fvl(loc,b)
  k = k + 1

  if ( la < lb ) then

    da = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + &
      nrmlc(3)*vcl(3,la) - nrmlc(4)

    if ( abs(da) <= dtol ) then
      ca = 0
    else if ( abs(da) <= mindis ) then
      return
    else if ( da < 0.0D+00 ) then
      ca = -1
    else
      ca = 1
    end if

    db = nrmlc(1)*vcl(1,lb) + nrmlc(2)*vcl(2,lb) + &
      nrmlc(3)*vcl(3,lb) - nrmlc(4)

    if ( abs(db) <= dtol ) then
      cb = 0
    else if ( abs(db) <= mindis ) then
      return
    else if ( db < 0.0D+00 ) then
      cb = -1
    else
      cb = 1
    end if

    if ( ca * cb <= 0 ) then

      nedg = nedg + 1
      ne = ne + 3

      if ( maxiw < ne ) then
        ierr = 6
        return
      end if

      if ( ccw == succ ) then
        iwk(ne-2) = a
        iwk(ne-1) = fvl(edgc,a)
      else
        iwk(ne-2) = b
        iwk(ne-1) = fvl(edgc,b)
      end if

      iwk(ne) = (ca+2)*10 + (cb+2)

    end if

  end if

  a = b

  if ( a /= facep(1,f) ) then
    go to 40
  end if

  nv = max ( nv, k )
  ptr = pfl(2,ptr)

  if ( ptr /= hfl(p) ) then
    go to 30
  end if

  if ( nedg <= 2 ) then
    ierr = 347
    return
  else if ( maxiw < ne + max ( nv+1, nedg+nedg+2 ) ) then
    ierr = 6
    return
  else if ( maxwk < nedg ) then
    ierr = 7
    return
  else if ( maxvc < nvc + 2 ) then
    ierr = 14
    return
  end if
!
!  Determine first intersection point of ray PT + T*DIR, 0 < T,
!  with boundary of polyhedron P.
!
  nep1 = ne + 1
  fmin = 0
  tmin = 0.0D+00
  ptr = hfl(p)

50 continue

  sf = pfl(1,ptr)
  f = abs ( sf )
  nflag = ( abs(facep(2,f) ) /= p )
  dotp = dot_product ( dir(1:3), nrml(1:3,f) )

  if ( nflag ) then
    dotp = -dotp
  end if

  if ( dotp <= tol ) then
    go to 90
  end if

  j = facep(1,f)
  la = fvl(loc,j)

  da = dot_product ( nrml(1:3,f), vcl(1:3,la) )

  t = ( da - dot_product ( nrml(1:3,f), pt(1:3) ) ) /dotp

  if ( nflag ) then
    t = -t
  end if

  if ( t <= tol .or. fmin /= 0 .and. tmin <= t ) then
    go to 90
  end if

  ipt(1:3) = pt(1:3) + t*dir(1:3)

  if ( sf < 0 .eqv. nflag ) then
    link = succ
  else
    link = pred
  end if

  i = nep1
  iwk(i) = la

  do

    j = fvl(link,j)
    i = i + 1
    iwk(i) = fvl(loc,j)
    if ( j == facep(1,f) ) then
      exit
    end if

  end do

  call ptpolg(3,3,i-nep1,1,iwk(nep1),vcl,ipt,nrml(1,f),dtol,inout)

  if ( 0 <= inout ) then
    fmin = f
    iomin = inout
    tmin = t
  end if

90 continue

  ptr = pfl(2,ptr)

  if ( ptr /= hfl(p) ) then
    go to 50
  end if

  if ( fmin == 0 ) then
    ierr = 347
    return
  end if

  ipt(1:3) = pt(1:3) + tmin * dir(1:3)

  if ( abs ( facep(2,fmin) ) == p ) then
    sf = facep(2,fmin)
  else
    sf = facep(3,fmin)
  end if

  if ( iomin == 1 ) then
    ibeg = 0
    iend = 1
    vbeg = facep(1,fmin)
    vend = vbeg
    go to 220
  end if
!
!  If IPT lies on boundary of face, determine edge containing IPT
!  and whether edge lies on cut plane.
!
  a = facep(1,fmin)

120 continue

  la = fvl(loc,a)

  do i = 1, 3
    cmax = max ( abs ( vcl(i,la) ), abs ( ipt(i) ) )
    if ( tol * cmax < abs ( vcl(i,la) - ipt(i) ) .and. tol < cmax ) then
      go to 140
    end if
  end do

  ierr = 348
  return

140 continue

  a = fvl(succ,a)
  if ( a /= facep(1,fmin) ) then
    go to 120
  end if

  cmin = 2.0D+00
  a = facep(1,fmin)
  la = fvl(loc,a)
  va(1:3) = vcl(1:3,la) - ipt(1:3)
  da = va(1)**2 + va(2)**2 + va(3)**2

160 continue

  b = fvl(succ,a)
  lb = fvl(loc,b)
  vb(1:3) = vcl(1:3,lb) - ipt(1:3)
  db = vb(1)**2 + vb(2)**2 + vb(3)**2
  dotp = dot_product ( va(1:3), vb(1:3) ) / sqrt ( da * db )

  if ( dotp <= -1.0D+00 + tol ) then
    j = a
    go to 190
  else if ( dotp < cmin ) then
    cmin = dotp
    j = a
  end if

  a = b
  la = lb
  da = db
  va(1:3) = vb(1:3)

  if ( a /= facep(1,fmin) ) then
    go to 160
  end if

190 continue

  la = fvl(loc,j)
  lb = fvl(loc,fvl(succ,j))
  da = dot_product ( nrmlc(1:3), vcl(1:3,la) ) - nrmlc(4)
!
!  Starting edge lies on cut plane.  Initialize fields for CUTFAC.
!
  if ( abs(da) <= dtol ) then

    k = 1
    if ( abs(nrmlc(1)) < abs(nrmlc(2)) ) then
      k = 2
    end if

    if ( abs(nrmlc(k)) < abs(nrmlc(3)) ) then
      k = 3
    end if

    if ( k == 1 ) then
      cmax = dir(2)*(vcl(3,lb) - vcl(3,la)) - dir(3)* (vcl(2,lb) - vcl(2,la))
    else if ( k == 2 ) then
      cmax = dir(3)*(vcl(1,lb) - vcl(1,la)) - dir(1)* (vcl(3,lb) - vcl(3,la))
    else
      cmax = dir(1)*(vcl(2,lb) - vcl(2,la)) - dir(2)* (vcl(1,lb) - vcl(1,la))
    end if

    if ( 0.0D+00 < cmax * nrmlc(k) ) then

      iwk(nep1) = la
      iwk(nep1+2) = lb

      if ( 0 < sf ) then

        if ( 0 < lb - la ) then
          ang = eang(j)
          j = fvl(edgc,j)
        else
          j = fvl(edga,j)
          ang = eang(j)
        end if

      else

        if ( lb - la < 0 ) then
          ang = eang(j)
        else
          ang = eang(fvl(edga,j))
        end if

      end if

    else

      iwk(nep1) = lb
      iwk(nep1+2) = la

      if ( sf < 0 ) then

        if ( lb - la < 0 ) then
          ang = eang(j)
          j = fvl(edgc,j)
        else
          j = fvl(edga,j)
          ang = eang(j)
        end if

      else

        if ( 0 < lb - la ) then
          ang = eang(j)
        else
          ang = eang(fvl(edga,j))
        end if

      end if

    end if

    iwk(nep1+3) = j
    f = fvl(facn,j)
    dotp = dot_product ( nrmlc(1:3), nrml(1:3,f) )

    if ( 1.0D+00 - tol < abs ( dotp ) ) then
      dotp = sign ( 1.0D+00, dotp )
    end if

    nflag = ( abs(facep(2,f)) /= p )

    if ( nflag ) then
      dotp = -dotp
    end if

    angn = pi - acos ( dotp )
    dotp = dot_product ( dir(1:3), nrml(1:3,f) )

    if ( nflag ) then
      dotp = -dotp
    end if

    if ( dotp < 0.0D+00 ) then
      angn = 2.0D+00 * pi - angn
    end if

    wk(1) = ang - angn

    if ( angn < angacc .or. wk(1) < angacc ) then
      return
    end if

    i = 1

200 continue

    if ( iwk(i) /= j .and. iwk(i+1) /= j ) then
      i = i + 3
      go to 200
    end if

    iwk(i) = iwk(ne-2)
    iwk(i+1) = iwk(ne-1)
    iwk(i+2) = iwk(ne)
    nedg = nedg - 1
    go to 300

  end if
!
!  Edge containing boundary point IPT does not lie on cut plane.
!
  dotp = dot_product (  nrmlc(1:3), vcl(1:3,lb) - vcl(1:3,la) )

  if ( sf < 0 ) then
    dotp = -dotp
  end if

  if ( 0.0D+00 < dotp ) then
    ibeg = 0
    iend = 0
    iwk(nep1+3) = -j
  else
    ibeg = 1
    iend = 1
    iwk(nep1+1) = j
  end if

  k = 1 - ibeg
  l(k) = nvc + k + 1
  vcl(1:3,l(k)) = ipt(1:3)
  vbeg = fvl(succ,j)
  vend = j
!
!  Find 1 or 2 points of intersection of face FMIN with cut plane.
!
220 continue

  ray(1) = nrmlc(2)*nrml(3,fmin) - nrmlc(3)*nrml(2,fmin)
  ray(2) = nrmlc(3)*nrml(1,fmin) - nrmlc(1)*nrml(3,fmin)
  ray(3) = nrmlc(1)*nrml(2,fmin) - nrmlc(2)*nrml(1,fmin)

  if ( abs(facep(2,fmin)) /= p ) then
    ray(1:3) = -ray(1:3)
  end if

  n = nvc + ibeg + 1
  ineg = 0
  ipos = 0
  smin = 0.0D+00
  k = 1
  if ( abs(ray(1)) < abs(ray(2)) ) k = 2
  if ( abs(ray(k)) < abs(ray(3)) ) k = 3
  a = vbeg
  la = fvl(loc,a)
  da = nrmlc(1)*vcl(1,la) + nrmlc(2)*vcl(2,la) + &
    nrmlc(3)*vcl(3,la) - nrmlc(4)

230 continue

  b = fvl(succ,a)
  lb = fvl(loc,b)
  db = nrmlc(1)*vcl(1,lb) + nrmlc(2)*vcl(2,lb) + &
    nrmlc(3)*vcl(3,lb) - nrmlc(4)

  if ( abs(da) <= dtol ) then

    t = (vcl(k,la) - ipt(k)) / ray(k)

    if ( t < 0.0D+00 ) then

      if ( ibeg == 0 ) then

        if ( ineg == 0 .or. smin < t ) then
          ineg = a
          smin = t
          l(0) = la
        end if

      end if

    else

      if ( iend == 1 ) then

        if ( ipos == 0 .or. t < tmin ) then
          ipos = a
          tmin = t
          l(1) = la
        end if

      end if
    end if

  else if ( abs ( db ) <= dtol ) then

    go to 250

  else if ( da * db < 0.0D+00 ) then

    va(1:3) = vcl(1:3,la) - ipt(1:3)
    vb(1:3) = vcl(1:3,la) - vcl(1:3,lb)
    cp(1) = ray(2)*vb(3) - ray(3)*vb(2)
    cp(2) = ray(3)*vb(1) - ray(1)*vb(3)
    cp(3) = ray(1)*vb(2) - ray(2)*vb(1)
    j = 1
    if ( abs(cp(1)) < abs(cp(2))) then
      j = 2
    end if

    if ( abs(cp(j)) < abs(cp(3))) then
      j = 3
    end if

    if ( j == 1 ) then
      t = (va(2)*vb(3) - va(3)*vb(2))/cp(1)
    else if ( j == 2 ) then
      t = (va(3)*vb(1) - va(1)*vb(3))/cp(2)
    else
      t = (va(1)*vb(2) - va(2)*vb(1))/cp(3)
    end if

    if ( t < 0.0D+00 ) then

      if ( ibeg == 0 ) then

        if ( ineg == 0 .or. smin < t ) then
          ineg = a
          smin = t
          l(0) = n
        end if

      end if

    else

      if ( iend == 1 ) then

        if ( ipos == 0 .or. t < tmin ) then
          ipos = a
          tmin = t
          l(1) = n
        end if

      end if

    end if

  end if

250 continue

  a = b
  la = lb
  da = db
  if ( a /= vend ) go to 230

  if ( ( ibeg == 0 .and. ineg == 0 ) .or. &
       ( iend == 1 .and. ipos == 0 ) ) then
    ierr = 349
    return
  end if

  if ( ibeg == 0 ) then
    if ( nvc < l(0) ) then
      iwk(nep1+1) = ineg
      vcl(1:3,n) = ipt(1:3) + smin * ray(1:3)
      n = n + 1
    end if
  end if

  if ( iend == 1 ) then

    if ( nvc < l(1) ) then
      l(1) = n
      iwk(nep1+3) = -ipos
      vcl(1:3,n) = ipt(1:3) + tmin*ray(1:3)
    else if ( sf < 0 ) then
      iwk(nep1+3) = -ipos
    else
      iwk(nep1+3) = -fvl(pred,ipos)
    end if

  end if
!
!  Set or update fields for routine CUTFAC.
!
  iwk(nep1) = l(0)
  iwk(nep1+2) = l(1)

  if ( iend == 0 .and. l(0) <= nvc ) then

    iwk(nep1+2) = nvc + 1
    vcl(1:3,nvc+1) = vcl(1:3,nvc+2)

  else if ( nvc < l(0) ) then

    a = iwk(nep1+1)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))

    if ( 0 < (lb - la) * sf ) then
      iwk(nep1+1) = fvl(edgc,a)
    else
      iwk(nep1+1) = fvl(edga,a)
    end if

  end if

  if ( nvc < l(1) ) then

    j = -iwk(nep1+3)
    i = 1

290 continue

    if ( iwk(i) /= j .and. iwk(i+1) /= j ) then
      i = i + 3
      go to 290
    end if

    iwk(i) = iwk(ne-2)
    iwk(i+1) = iwk(ne-1)
    iwk(i+2) = iwk(ne)
    nedg = nedg - 1

  end if

  dotp = dot_product ( nrmlc(1:3), nrml(1:3,fmin) )

  if ( 1.0D+00 - tol < abs ( dotp ) ) then
    dotp = sign ( 1.0D+00, dotp )
  end if

  if ( abs(facep(2,fmin)) /= p ) then
    dotp = -dotp
  end if

  wk(1) = acos(dotp)
  angn = pi - wk(1)

  if ( angn < angacc .or. wk(1) < angacc ) then
    return
  end if

300 continue

  call cutfac(p,nrmlc,angacc,dtol,nvc,maxvc,vcl,facep,nrml,fvl,eang, &
     nedg,iwk,nce,iwk(nep1),wk,rflag,ierr)

  if ( ierr /= 0 ) then
    return
  end if

  if ( rflag ) then
    call insfac(p,nrmlc,nce,iwk(nep1),wk,nvc,nface,nvert,npolh, &
      npf,maxfp,maxfv,maxhf,maxpf,vcl,facep,factyp,nrml,fvl,eang, &
      hfl,pfl,ierr)
  end if

  return
end
