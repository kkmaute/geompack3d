subroutine sepfac ( p, cntr, nrmlc, angacc, mxcos, dtol, nvc, maxvc, vcl, &
  facep, nrml, fvl, eang, nce, cedge, cdang, aflag, ierr )

!*****************************************************************************80
!
!! SEPFAC traces out a separator in a convex polyhedron.
!
!  Discussion:
!
!    This routine traces out a separator or cut face in a convex polyhedron P
!    given a starting edge on boundary and normal to cut face.
!
!    Accept a cut face if it creates no small angles or short edges.
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
!    Input, real ( kind = 8 ) CNTR(1:3), the weighted centroid of polyhedron.
!
!    Input, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut plane
!    plus right hand side constant term of plane equation.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) MXCOS, the maximum cosine of angle allowed
!    for angles subtended by new subedges with respect to centroid.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine
!    if point is on cut plane.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC+2), the vertex
!    coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:*), the edge angles.
!
!    Input, integer ( kind = 4 ) CEDGE(1:2,0:1); CEDGE(1,0) = LU, CEDGE(1,1) = LV where
!    LU, LV are indices of VCL and are vertices of starting edge;
!    LU may be NVC+1 and LV may be NVC+1 or NVC+2 to indicate
!    a point on interior of an edge; CEDGE(2,1) is specified
!    as described for output below; if NVC < LU, CEDGE(2,0)
!    takes on the value ABS(CEDGE(2,NCE)) described for output
!    below, else CEDGE(2,0) is not used.
!
!    Input, real ( kind = 8 ) CDANG(1), dihedral angle at starting edge
!    determined by cut plane in positive half-space.
!
!    Output, integer ( kind = 4 ) NCE, the number of edges in cut polygon (if acceptable);
!    it is assumed there is enough space in the following two arrays.
!
!    Output, integer ( kind = 4 ) CEDGE(1:2,0:NCE); CEDGE(1,I) is an index of VCL,
!    NVC < indices are new points; CEDGE(2,I) = J indicates that edge of cut
!    face ending at CEDGE(1,I) is edge from J to FVL(SUCC,J)
!    if 0 < J; else if J < 0 then edge of cut face ending at
!    CEDGE(1,I) is a new edge and CEDGE(1,I) lies on edge from
!    -J to FVL(SUC,-J) and new edge lies in face FVL(FACN,-J);
!    CEDGE(2,I) always refers to an edge in the subpolyhedron
!    in negative half-space; CEDGE(1,NCE) = CEDGE(1,0);
!    CEDGE(2,0) is not used.
!
!    Output, real ( kind = 8 ) CDANG(1:NCE), the dihedral angles created by
!    edges of cut polygon in positive half-space; negative sign for angle
!    I indicates that face containing edge I is oriented CW in polyhedron P.
!
!    Output, logical AFLAG, TRUE iff separator face is acceptable.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxvc

  integer ( kind = 4 ) a
  logical              aflag
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angr
  integer ( kind = 4 ) b
  integer ( kind = 4 ) ccwfl
  real    ( kind = 8 ) cdang(*)
  integer ( kind = 4 ) cedge(2,0:*)
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) cp
  real    ( kind = 8 ) de(3)
  real    ( kind = 8 ) dee(3)
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dist
  real    ( kind = 8 ) distb
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dtol
  integer ( kind = 4 ) e
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) ee
  logical              eflag
  integer ( kind = 4 ) estrt
  integer ( kind = 4 ) estop
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fl
  integer ( kind = 4 ) fr
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lu
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) lw1
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mxcos
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  real    ( kind = 8 ) nmax
  real    ( kind = 8 ) nrml(3,*)
  real    ( kind = 8 ) nrmlc(4)
  real    ( kind = 8 ) ntol
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) sgn
  integer ( kind = 4 ) sp
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) w

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  aflag = .false.
  n = max ( nvc, cedge(1,0), cedge(1,1) )
  nce = 1
  lu = cedge(1,0)
  lw = lu
  lw1 = cedge(1,1)
  w = cedge(2,1)

  if ( 0 < w ) then
    fl = fvl(facn,w)
    if ( lw1 < lw ) then
      fr = fvl(facn,fvl(edgc,w))
    else
      fr = fvl(facn,fvl(edga,w))
    end if
  else
    fl = fvl(facn,-w)
    fr = fl
  end if

  dir(1:3) = vcl(1:3,lw1) - vcl(1:3,lw)

  if ( abs(facep(2,fl)) == p ) then
    ccwfl = facep(2,fl)
  else
    ccwfl = facep(3,fl)
  end if

  if ( ccwfl < 0 ) then
    cdang(1) = -cdang(1)
  end if
!
!  LW, LW1, FL, FR, DIR are updated before each iteration of loop.
!  counterclockwiseFL = P (-P) if FL is CCW (CW) according to SUCC traversal.
!
10 continue
!
!  LW1 is new vertex on interior of edge E. FL is used for
!  previous and next faces containing edges of cut polygon.
!
  if ( nvc < lw1 ) then

    e = -cedge(2,nce)
    la = fvl(loc,e)
    lb = fvl(loc,fvl(succ,e))

    if ( 0 < (lb - la)*ccwfl ) then
      ee = fvl(edgc,e)
    else
      ee = fvl(edga,e)
    end if

    fl = fvl(facn,ee)

    if ( abs(facep(2,fl)) == p ) then
      ccwfl = facep(2,fl)
    else
      ccwfl = facep(3,fl)
    end if

    dir(1) = nrmlc(2)*nrml(3,fl) - nrmlc(3)*nrml(2,fl)
    dir(2) = nrmlc(3)*nrml(1,fl) - nrmlc(1)*nrml(3,fl)
    dir(3) = nrmlc(1)*nrml(2,fl) - nrmlc(2)*nrml(1,fl)

    if ( abs(facep(2,fl)) /= p ) then
      dir(1:3) = -dir(1:3)
    end if

  else
!
!  LW1 is existing vertex of polyhedron P. FL (FL and FR) is
!  previous face if edge ending at LW1 is new (already exists).
!  In former case, -CEDGE(2,NCE) is the edge of FL incident on
!  LW1 which will lie only in subpolyhedron PL. In latter case,
!  CEDGE(2,NCE) is an edge of FL. Cycle thru edges, faces counterclockwise
!  (from outside) between edges ESTRT, ESTOP.
!
    eflag = ( 0 < cedge(2,nce) )

    if ( .not. eflag ) then

      estrt = -cedge(2,nce)
      sp = ccwfl

      if ( 0 < ccwfl ) then
        estop = fvl(succ,estrt)
      else
        estop = fvl(pred,estrt)
      end if

    else

      w = cedge(2,nce)

      if ( 0 < ccwfl ) then
        estrt = fvl(pred,w)
        l = lw - lw1
      else
        estrt = fvl(succ,w)
        l = lw1 - lw
      end if

      if ( 0 < l*ccwfl ) then
        w = fvl(edgc,w)
      else
        w = fvl(edga,w)
      end if

      if ( abs(facep(2,fr)) == p ) then
        sp = facep(2,fr)
      else
        sp = facep(3,fr)
      end if

      if ( 0 < sp ) then
        estop = fvl(succ,w)
      else
        estop = fvl(pred,w)
      end if

    end if

    la = fvl(loc,estop)
    lb = fvl(loc,fvl(succ,estop))

    if ( 0 < (lb - la)*sp ) then
      estop = fvl(edgc,estop)
    else
      estop = fvl(edga,estop)
    end if

    e = estrt

20  continue

    if ( eflag .or. ( e /= estrt .and. e /= estop ) ) then
!
!  Determine if edge lies on cut plane.
!
      if ( fvl(loc,e) == lw1 ) then
        l = fvl(loc,fvl(succ,e))
      else
        l = fvl(loc,e)
      end if

      dist = nrmlc(1)*vcl(1,l) + nrmlc(2)*vcl(2,l) + &
        nrmlc(3)*vcl(3,l) - nrmlc(4)

      if ( abs(dist) <= dtol ) then

        dir(1:3) = vcl(1:3,l) - vcl(1:3,lw1)
        lw = lw1
        lw1 = l
        nce = nce + 1
        cedge(1,nce) = lw1
        cedge(2,nce) = e
        fl = fvl(facn,e)

        if ( abs(facep(2,fl)) == p ) then
          ccwfl = facep(2,fl)
        else
          ccwfl = facep(3,fl)
        end if

        if ( 0 < ccwfl ) then
          l = lw - lw1
        else
          l = lw1 - lw
        end if

        if ( 0 < l*ccwfl ) then
          fr = fvl(facn,fvl(edgc,e))
        else
          fr = fvl(facn,fvl(edga,e))
        end if

        go to 40

      end if

    end if

    if ( e == estop ) then
      ierr = 334
      return
    end if
!
!  Determine if cut plane intersects interior of face.
!
    la = fvl(loc,e)
    lb = fvl(loc,fvl(succ,e))

    if ( 0 < (lb - la)*ccwfl ) then
      e = fvl(edgc,e)
    else
      e = fvl(edga,e)
    end if

    fl = fvl(facn,e)

    if ( abs(facep(2,fl)) == p ) then
      ccwfl = facep(2,fl)
    else
      ccwfl = facep(3,fl)
    end if

    if ( 0 < ccwfl ) then

      ee = fvl(pred,e)
      la = fvl(loc,fvl(succ,e))
      lb = fvl(loc,ee)

    else

      ee = fvl(succ,e)
      la = fvl(loc,e)
      lb = fvl(loc,fvl(succ,ee))

    end if

    dir(1) = nrmlc(2)*nrml(3,fl) - nrmlc(3)*nrml(2,fl)
    dir(2) = nrmlc(3)*nrml(1,fl) - nrmlc(1)*nrml(3,fl)
    dir(3) = nrmlc(1)*nrml(2,fl) - nrmlc(2)*nrml(1,fl)
    sgn = 1

    if ( abs(facep(2,fl)) /= p ) then
      dir(1:3) = -dir(1:3)
      sgn = -1
    end if

    k = 1
    if ( abs(nrml(1,fl)) < abs(nrml(2,fl)) ) then
      k = 2
    end if

    if ( abs(nrml(k,fl)) < abs(nrml(3,fl)) ) then
      k = 3
    end if

    nmax = sgn * nrml(k,fl)

    de(1:3) = vcl(1:3,la) - vcl(1:3,lw1)
    dee(1:3) = vcl(1:3,lb) - vcl(1:3,lw1)

    ntol = tol * max ( abs ( de(1) ), abs ( de(2) ), abs( de(3) ), &
      abs ( dee(1) ), abs ( dee(2) ), abs ( dee(3) ) )
    e = ee

    if ( k == 1 ) then
      cp = de(2)*dir(3) - de(3)*dir(2)
    else if ( k == 2 ) then
      cp = de(3)*dir(1) - de(1)*dir(3)
    else
      cp = de(1)*dir(2) - de(2)*dir(1)
    end if

    if ( abs(cp) <= ntol .or. cp*nmax < 0.0D+00) go to 20

    if ( k == 1 ) then
      cp = dir(2)*dee(3) - dir(3)*dee(2)
    else if ( k == 2 ) then
      cp = dir(3)*dee(1) - dir(1)*dee(3)
    else
      cp = dir(1)*dee(2) - dir(2)*dee(1)
    end if

    if ( abs(cp) <= ntol .or. cp*nmax < 0.0D+00) go to 20

  end if
!
!  Next cut edge is in interior of face FL. Determine LW1.
!  Edge EE containing LW has been saved above in both cases.
!
  fr = fl
  lw = lw1

  if ( nvc < lw ) then
    a = fvl(succ,ee)
    e = ee
  else if ( fvl(loc,ee) == lw ) then
    a = fvl(succ,ee)
    e = fvl(pred,ee)
  else
    a = fvl(succ,fvl(succ,ee))
    e = ee
  end if

  la = fvl(loc,a)
  dist = nrmlc(4) - nrmlc(1)*vcl(1,la) - nrmlc(2)*vcl(2,la) - &
    nrmlc(3)*vcl(3,la)

30 continue

  b = fvl(succ,a)
  lb = fvl(loc,b)
  distb = nrmlc(4) - nrmlc(1)*vcl(1,lb) - nrmlc(2)*vcl(2,lb) - &
    nrmlc(3)*vcl(3,lb)

  if ( b /= e .and. abs ( distb ) <= dtol ) then

    lw1 = lb
    nce = nce + 1
    cedge(1,nce) = lw1

    if ( 0 < ccwfl ) then
      cedge(2,nce) = -a
    else
      cedge(2,nce) = -b
    end if

    go to 40

  else if ( dist*distb < 0.0D+00 ) then

    if ( nvc < lu ) then

      if ( a == cedge(2,0) ) then
        lw1 = lu
        nce = nce + 1
        cedge(1,nce) = lu
        cedge(2,nce) = -a
        go to 40
      end if

    end if

    de(1:3) = vcl(1:3,lb) - vcl(1:3,la)

    t = dist/(nrmlc(1)*de(1)+nrmlc(2)*de(2)+nrmlc(3)*de(3))
    n = n + 1

    if ( maxvc < n ) then
      ierr = 14
      return
    end if

    vcl(1:3,n) = vcl(1:3,la) + t * de(1:3)
!
!  If angles subtended from centroid by subedges are too
!  small, reject cut plane.
!
    de(1:3) = vcl(1:3,n) - cntr(1:3)
    leng = de(1)**2 + de(2)**2 + de(3)**2

    dee(1:3) = vcl(1:3,la) - cntr(1:3)

    t = dee(1)**2 + dee(2)**2 + dee(3)**2
    dotp = (de(1)*dee(1) + de(2)*dee(2) + de(3)*dee(3))/ sqrt(leng*t)

    if ( mxcos < dotp ) then

      if ( msglvl == 4 ) then
        write ( *,600) 'rejected due to short subedge'
      end if

      return

    end if

    dee(1:3) = vcl(1:3,lb) - cntr(1:3)

    t = dee(1)**2 + dee(2)**2 + dee(3)**2
    dotp = (de(1)*dee(1) + de(2)*dee(2) + de(3)*dee(3))/ sqrt(leng*t)

    if ( mxcos < dotp ) then

      if ( msglvl == 4 ) then
        write ( *, '(a)' ) '  Rejected due to short subedge'
      end if

      return

    end if

    lw1 = n
    nce = nce + 1
    cedge(1,nce) = n
    cedge(2,nce) = -a
    go to 40

  end if

  if ( b /= e ) then
    a = b
    la = lb
    dist = distb
    go to 30
  else
    ierr = 335
    return
  end if
!
!  Compute dihedral angles due to cut plane at cut edge. If any
!  angle is too small, reject cut plane.
!
40 continue

  dotp = dot_product ( nrmlc(1:3), nrml(1:3,fr) )

  if ( 1.0D+00 - tol < abs ( dotp ) ) then
    dotp = sign ( 1.0D+00, dotp )
  end if

  if ( abs(facep(2,fr)) /= p ) then
    dotp = -dotp
  end if

  angr = acos(dotp)

  if ( fr == fl ) then

    ang = pi

  else

    a = cedge(2,nce)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))

    if ( 0 < (lb - la)*ccwfl ) then
      ang = eang(a)
    else
      ang = eang(fvl(edga,a))
    end if

  end if

  if ( angr < angacc .or. ang - angr < angacc ) then
    if ( msglvl == 4 ) then
      write ( *,600) 'rejected due to small angle'
    end if
    return
  end if

  if ( 0 < ccwfl ) then
    cdang(nce) = angr
  else
    cdang(nce) = -angr
  end if

  if ( lw1 /= lu ) go to 10

  aflag = .true.

  if ( msglvl == 4 ) then
    write ( *,600) 'cedge(1:2), cdang'
    do k = 1,nce
      write ( *,610) k,cedge(1,k),cedge(2,k),cdang(k)*180.0D+00 / pi
    end do
  end if

  600 format (4x,a)
  610 format (4x,3i7,f12.5)

  return
end
