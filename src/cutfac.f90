subroutine cutfac ( p, nrmlc, angacc, dtol, nvc, maxvc, vcl, facep, nrml, &
  fvl, eang, nedgc, pedge, nce, cedge, cdang, rflag, ierr )

!*****************************************************************************80
!
!! CUTFAC traces a cut face of a polyhedron from a starting edge.
!
!  Discussion:
!
!    This routine traces out the cut face in polyhedron P given a starting
!    edge, for example, a reflex edge or an edge in the interior of a face.
!
!    It accepts the cut face if it creates no small angles, it is an outer
!    boundary and there are no holes (inner polygons) in it.  It is
!    assumed the starting edge does not lie on a double-occurring face.
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
!    Input, integer ( kind = 4 ) P, the polyhedron index.
!
!    Input, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut
!    plane plus right hand side constant term of plane equation.
!
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine
!    if point is on cut plane.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:*), the edge angles.
!
!    Input/output, integer ( kind = 4 ) NEDGC, the number of edges which intersect cut plane.
!
!    Input/output, integer ( kind = 4 ) PEDGE(1:3,1:NEDGC), the edges of polyhedron P that
!    intersect cut plane excluding starting edge if it is an edge of P or
!    the edge containing CEDGE(1,1) if NVC < CEDGE(1,1);
!    PEDGE(1,I), PEDGE(2,I) are indices of FVL; if PEDGE(1,I)
!    = A and B = FVL(SUCC,A ) then PEDGE(3,I) = 10*CA+CB where
!    CA = 1, 2, or 3 if vertex A in negative half-space, on
!    cut plane, or in positive half-space determined by cut
!    plane, and similarly for CB.
!
!    Input, integer ( kind = 4 ) CEDGE(1:2,0:1); CEDGE(1,0) = LV, CEDGE(1,1) = LU where
!    LU, LV are indices of VCL and are vertices of starting edge; if
!    called by RESEDG, LU < LV <= NVC are on reflex edge and
!    CEDGE(2,1) = U is index of FVL indicating reflex edge;
!    LV may be NVC+1 and LU may be NVC+1 or NVC+2 to indicate
!    a point on interior of an edge; CEDGE(2,1) is specified
!    as described for output below; if NVC < LV, CEDGE(2,0)
!    takes on the value ABS ( CEDGE(2,NCE) ) described for output
!    below, else CEDGE(2,0) is not used.
!
!    Input, real ( kind = 8 ) CDANG(1), the dihedral angle at starting edge
!    determined by cut plane in positive half-space.
!
!    Output, integer ( kind = 4 ) NCE, the number of edges in cut polygon; it is assumed
!    there is enough space in the following two arrays.
!
!    Output, integer ( kind = 4 ) CEDGE(1:2,0:NCE); CEDGE(1,I) is an index of VCL,
!    points with indices greater than NVC are new points; CEDGE(2,I) = J
!    indicates that edge of cut face ending at CEDGE(1,I) is edge from J
!    to FVL(SUCC,J) if J is greater than 0; else if J < 0 then edge of cut
!    face ending at CEDGE(1,I) is a new edge and CEDGE(1,I) lies on edge from
!    -J to FVL(SUC,-J) and new edge lies in face FVL(FACN,-J);
!    CEDGE(2,I) always refers to an edge in the subpolyhedron
!    in negative half-space; CEDGE(1,NCE) = CEDGE(1,0).
!
!    Output, real ( kind = 8 ) CDANG(1:NCE), the dihedral angles created by
!    edges of cut polygon in positive half-space; negative sign for angle I
!    indicates that face containing edge I is oriented clockwise in
!    polyhedron P.
!
!    Output, logical RFLAG, TRUE iff reflex or starting edge is resolved.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ), parameter :: maxev = 20
  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) nedgc

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) angr
  integer ( kind = 4 ) ca
  integer ( kind = 4 ) cb
  integer ( kind = 4 ) ccwfl
  real    ( kind = 8 ) cdang(*)
  integer ( kind = 4 ) cedge(2,0:*)
  real    ( kind = 8 ) cmax
  real    ( kind = 8 ) cp(3)
  real    ( kind = 8 ) de(3)
  real    ( kind = 8 ) dee(3)
  real    ( kind = 8 ) dir(3)
  real    ( kind = 8 ) dir1(3)
  real    ( kind = 8 ) dirsq
  real    ( kind = 8 ) dir1sq
  real    ( kind = 8 ) dist
  logical              dof
  real    ( kind = 8 ) dotp
  real    ( kind = 8 ) dsave(3)
  real    ( kind = 8 ) dtol
  integer ( kind = 4 ) e
  real    ( kind = 8 ) eang(*)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) ee
  logical              eflag
  integer ( kind = 4 ) estrt
  integer ( kind = 4 ) estop
  integer ( kind = 4 ) ev(maxev)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fl
  integer ( kind = 4 ) fr
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iamin
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) imin
  integer ( kind = 4 ) inout
  real    ( kind = 8 ) intang
  integer ( kind = 4 ) isave
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kmax
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lw
  integer ( kind = 4 ) lw1
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  integer ( kind = 4 ) nev
  real    ( kind = 8 ) nmax
  real    ( kind = 8 ) nrml(3,*)
  real    ( kind = 8 ) nrmlc(4)
  real    ( kind = 8 ) ntol
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  integer ( kind = 4 ) pedge(3,nedgc)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pi2
  integer ( kind = 4 ), parameter :: pred = 4
  logical              rflag
  real    ( kind = 8 ) rhs(3)
  real    ( kind = 8 ) s
  integer ( kind = 4 ) sf
  integer ( kind = 4 ) sgn
  integer ( kind = 4 ) sp
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tmin
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
  integer ( kind = 4 ) w

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  rflag = .false.
  pi2 = 2.0D+00 * pi
  n = max ( nvc, cedge(1,0), cedge(1,1) )
  nce = 1
  w = cedge(2,1)
  lv = cedge(1,0)
  lw = lv
  lw1 = cedge(1,1)

  if ( lw1 <= nvc ) then
    nev = 1
    ev(nev) = lw1
  else
    nev = 0
  end if

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

  dir(1) = vcl(1,lw1) - vcl(1,lw)
  dir(2) = vcl(2,lw1) - vcl(2,lw)
  dir(3) = vcl(3,lw1) - vcl(3,lw)

  kmax = 1
  if ( abs ( nrmlc(1) ) < abs ( nrmlc(2) ) ) then
    kmax = 2
  end if

  if ( abs ( nrmlc(kmax) ) < abs ( nrmlc(3) ) ) then
    kmax = 3
  end if

  if ( abs ( facep(2,fl) ) == p ) then
    ccwfl = facep(2,fl)
  else
    ccwfl = facep(3,fl)
  end if

  if ( ccwfl < 0 ) then
    cdang(1) = -cdang(1)
  end if
!
!  LW, LW1, FL, FR, DIR are updated before each iteration of loop.
!  CCWFL = P (-P) if FL is counterclockwise (clockwise) according
!  to SUCC traversal.
!
10 continue

  if ( nvc < lw1 ) then
!
!  LW1 is new vertex on interior of edge E.  FL is used for
!  previous and next faces containing edges of cut polygon.
!
    e = -cedge(2,nce)
    la = fvl(loc,e)
    lb = fvl(loc,fvl(succ,e))

    if ( 0 < ( lb - la ) * ccwfl ) then
      ee = fvl(edgc,e)
    else
      ee = fvl(edga,e)
    end if

    fl = fvl(facn,ee)
    dof = ( abs ( facep(2,fl) ) == abs ( facep(3,fl) ) )

    if ( dof ) then
      l = fvl(loc,ee)
      if ( l == la ) then
        ccwfl = -ccwfl
      end if
    else if ( abs ( facep(2,fl) ) == p ) then
      ccwfl = facep(2,fl)
    else
      ccwfl = facep(3,fl)
    end if

    dir(1) = nrmlc(2) * nrml(3,fl) - nrmlc(3) * nrml(2,fl)
    dir(2) = nrmlc(3) * nrml(1,fl) - nrmlc(1) * nrml(3,fl)
    dir(3) = nrmlc(1) * nrml(2,fl) - nrmlc(2) * nrml(1,fl)

    if ( abs ( facep(2,fl) ) /= p .or. ( dof .and. ccwfl < 0 ) ) then
      dir(1:3) = -dir(1:3)
    end if

    go to 70

  else
!
!  LW1 is existing vertex of polyhedron P. FL (FL and FR) is
!  previous face if edge ending at LW1 is new (already exists).
!  In former case, -CEDGE(2,NCE) is the edge of FL incident on
!  LW1 which will lie only in subpolyhedron PL.  In latter case,
!  CEDGE(2,NCE) is an edge of FL.  Cycle through edges, faces counterclockwise
!  (from outside) between edges ESTRT, ESTOP.
!  If LW1 lies on a doubly-occurring face, there are 2 cycles
!  around LW1 and the correct one is chosen based on CCWFL.
!
    iamin = pi2
    imin = 0
    dirsq = dir(1)**2 + dir(2)**2 + dir(3)**2
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
      la = fvl(loc,w)

      if ( 0 < ccwfl ) then
        estrt = fvl(pred,w)
      else
        estrt = fvl(succ,w)
      end if

      if ( la == lw ) then
        l = lw1 - lw
      else
        l = lw - lw1
      end if

      if ( 0 < l*ccwfl ) then
        w = fvl(edgc,w)
      else
        w = fvl(edga,w)
      end if

      if ( abs ( facep(2,fr) ) == abs ( facep(3,fr) ) ) then

        lb = fvl(loc,w)

        if ( la == lb ) then
          sp = -ccwfl
        else
          sp = ccwfl
        end if

      else if ( abs ( facep(2,fr) ) == p ) then

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
    sf = ccwfl

20  continue

    if ( eflag .or. ( e /= estrt .and. e /= estop ) ) then

      if ( fvl(loc,e) == lw1 ) then
        l = fvl(loc,fvl(succ,e))
      else
        l = fvl(loc,e)
      end if

      dist = nrmlc(1)*vcl(1,l) + nrmlc(2)*vcl(2,l) + &
        nrmlc(3)*vcl(3,l) - nrmlc(4)

      if ( abs ( dist ) <= dtol ) then

        dir1(1) = vcl(1,l) - vcl(1,lw1)
        dir1(2) = vcl(2,l) - vcl(2,lw1)
        dir1(3) = vcl(3,l) - vcl(3,lw1)
        dir1sq = dir1(1)**2 + dir1(2)**2 + dir1(3)**2
        dotp = -(dir(1)*dir1(1) + dir(2)*dir1(2) + dir(3)* &
          dir1(3))/sqrt(dirsq*dir1sq)

        if ( 1.0D+00 - tol < abs ( dotp ) ) then
          dotp = sign ( 1.0D+00, dotp )
        end if

        if ( kmax == 1 ) then
          cp(1) = dir(2)*dir1(3) - dir(3)*dir1(2)
        else if ( kmax == 2 ) then
          cp(2) = dir(3)*dir1(1) - dir(1)*dir1(3)
        else
          cp(3) = dir(1)*dir1(2) - dir(2)*dir1(1)
        end if

        if ( abs ( cp(kmax) ) <= tol * max ( abs ( dir(1) ), &
          abs ( dir(2) ), abs ( dir(3) ), abs ( dir1(1) ), &
          abs ( dir1(2) ), abs ( dir1(3) ) ) ) then

          intang = pi

        else if ( 0.0D+00 < cp(kmax) * nrmlc(kmax) ) then
          intang = acos(dotp)
        else
          intang = pi2 - acos ( dotp )
        end if

        if ( intang < iamin ) then
          iamin = intang
          imin = e
          dsave(1) = dir1(1)
          dsave(2) = dir1(2)
          dsave(3) = dir1(3)
        end if

      end if

    end if

    if ( e == estop ) then
      go to 40
    end if

    la = fvl(loc,e)
    lb = fvl(loc,fvl(succ,e))

    if ( 0 < (lb - la)*sf ) then
      e = fvl(edgc,e)
    else
      e = fvl(edga,e)
    end if

    f = fvl(facn,e)
    dof = ( abs ( facep(2,f) ) == abs ( facep(3,f) ) )

    if ( dof ) then
      l = fvl(loc,e)
      if ( l == la ) then
        sf = -sf
      end if
    else if ( abs ( facep(2,f) ) == p ) then
      sf = facep(2,f)
    else
      sf = facep(3,f)
    end if

    if ( 0 < sf ) then
      ee = fvl(pred,e)
      la = fvl(loc,fvl(succ,e))
      lb = fvl(loc,ee)
    else
      ee = fvl(succ,e)
      la = fvl(loc,e)
      lb = fvl(loc,fvl(succ,ee))
    end if

    dir1(1) = nrmlc(2) * nrml(3,f) - nrmlc(3) * nrml(2,f)
    dir1(2) = nrmlc(3) * nrml(1,f) - nrmlc(1) * nrml(3,f)
    dir1(3) = nrmlc(1) * nrml(2,f) - nrmlc(2) * nrml(1,f)

    if ( max ( abs ( dir1(1) ), abs ( dir1(2) ), abs ( dir1(3) ) ) <= tol ) then
      go to 30
    end if

    sgn = 1

    if ( abs ( facep(2,f ) ) /= p .or. ( dof .and. sf < 0 ) ) then
      dir1(1:3) = -dir1(1:3)
      sgn = -1
    end if

    k = 1

    if ( abs ( nrml(1,f) ) < abs ( nrml(2,f) ) ) then
      k = 2
    end if

    if ( abs ( nrml(k,f) ) < abs ( nrml(3,f) ) ) then
      k = 3
    end if

    nmax = sgn * nrml(k,f)

    de(1:3) = vcl(1:3,la) - vcl(1:3,lw1)
    dee(1:3) = vcl(1:3,lb) - vcl(1:3,lw1)

    ntol = tol * max ( abs ( de(1) ), abs ( de(2) ), abs ( de(3) ), &
      abs ( dee(1) ), abs ( dee(2) ), abs ( dee(3) ) )

    if ( k == 1 ) then
      cp(1) = de(2) * dee(3) - de(3) * dee(2)
    else if ( k == 2 ) then
      cp(2) = de(3) * dee(1) - de(1) * dee(3)
    else
      cp(3) = de(1) * dee(2) - de(2) * dee(1)
    end if

    if ( abs ( cp(k) ) <= ntol .or. 0.0D+00 < cp(k) * nmax ) then
      if ( k == 1) cp(1) = de(2)*dir1(3) - de(3)*dir1(2)
      if ( k == 2) cp(2) = de(3)*dir1(1) - de(1)*dir1(3)
      if ( k == 3) cp(3) = de(1)*dir1(2) - de(2)*dir1(1)

      if ( abs ( cp(k) ) <= ntol .or. cp(k) * nmax < 0.0D+00 ) then
        go to 30
      end if

      if ( k == 1 ) then
        cp(1) = dir1(2) * dee(3) - dir1(3) * dee(2)
      else if ( k == 2 ) then
        cp(2) = dir1(3) * dee(1) - dir1(1) * dee(3)
      else if ( k == 3 ) then
        cp(3) = dir1(1) * dee(2) - dir1(2) * dee(1)
      end if

      if ( abs ( cp(k) ) <= ntol .or. cp(k) * nmax < 0.0D+00 ) then
        go to 30
      end if

    else

      if ( k == 1 ) then
        cp(1) = dir1(2)*de(3) - dir1(3)*de(2)
      else if ( k == 2 ) then
        cp(2) = dir1(3)*de(1) - dir1(1)*de(3)
      else if ( k == 3 ) then
        cp(3) = dir1(1)*de(2) - dir1(2)*de(1)
      end if

      if ( abs ( cp(k) ) <= ntol .or. 0.0D+00 < cp(k) * nmax ) then

        if ( k == 1) cp(1) = dee(2)*dir1(3)-dee(3)*dir1(2)
        if ( k == 2) cp(2) = dee(3)*dir1(1)-dee(1)*dir1(3)
        if ( k == 3) cp(3) = dee(1)*dir1(2)-dee(2)*dir1(1)

        if ( abs(cp(k)) <= ntol .or. 0.0D+00 < cp(k)*nmax ) then
          go to 30
        end if

      end if

    end if

    dir1sq = dir1(1)**2 + dir1(2)**2 + dir1(3)**2
    dotp = - dot_product ( dir(1:3), dir1(1:3) ) / sqrt ( dirsq * dir1sq )

    if ( 1.0D+00 - tol < abs ( dotp ) ) then
      dotp = sign ( 1.0D+00, dotp )
    end if

    if ( kmax == 1 ) then
      cp(1) = dir(2)*dir1(3) - dir(3)*dir1(2)
    else if ( kmax == 2 ) then
      cp(2) = dir(3)*dir1(1) - dir(1)*dir1(3)
    else
      cp(3) = dir(1)*dir1(2) - dir(2)*dir1(1)
    end if

    if ( abs ( cp(kmax) ) <= tol * max ( abs ( dir(1) ), &
      abs ( dir(2) ), abs ( dir(3) ), abs ( dir1(1) ), &
      abs ( dir1(2) ), abs ( dir1(3) ) ) ) then

      intang = pi

    else if ( 0.0D+00 < cp(kmax) * nrmlc(kmax) ) then
      intang = acos(dotp)
    else
      intang = pi2 - acos(dotp)
    end if

    if ( intang < iamin ) then
      iamin = intang
      imin = -f
      ccwfl = sf
      dsave(1:3) = dir1(1:3)
    end if

30  continue

    e = ee
    go to 20

40  continue

    if ( imin == 0 ) then

      return

    else if ( 0 < imin ) then

      dir(1:3) = dsave(1:3)
      lw = lw1
      la = fvl(loc,imin)
      lb = fvl(loc,fvl(succ,imin))

      if ( la == lw1 ) then
        lw1 = lb
      else
        lw1 = la
      end if

      nce = nce + 1
      cedge(1,nce) = lw1
      cedge(2,nce) = imin
      fl = fvl(facn,imin)
      dof = ( abs(facep(2,fl) ) == abs ( facep(3,fl) ) )

      if ( dof ) then

        if ( la == lw ) then
          ccwfl = -p
        else
          ccwfl = p
        end if

      else if ( abs(facep(2,fl)) == p ) then

        ccwfl = facep(2,fl)

      else

        ccwfl = facep(3,fl)

      end if

      if ( 0 < (lb - la)*ccwfl ) then
        fr = fvl(facn,fvl(edgc,imin))
      else
        fr = fvl(facn,fvl(edga,imin))
      end if

      k = 1

50    continue

      if ( pedge(1,k) == imin .or. pedge(2,k) == imin ) then

        do i = 1, 3
          j = pedge(i,k)
          pedge(i,k) = pedge(i,nedgc)
          pedge(i,nedgc) = j
        end do

        nedgc = nedgc - 1

      else

        k = k + 1
        go to 50

      end if

      go to 110

    else

      dir(1) = dsave(1)
      dir(2) = dsave(2)
      dir(3) = dsave(3)
      fl = -imin
      go to 70

    end if

  end if
!
!  Determine LW1 from direction DIR in interior of face FL.
!
70 continue

  lw = lw1
  fr = 0
  imin = 0
  tmin = 0.0D+00
  k = 1

  if ( abs(dir(1)) < abs(dir(2)) ) then
    k = 2
  end if

  if ( abs(dir(k)) < abs(dir(3)) ) then
    k = 3
  end if

  ntol = tol * abs(dir(k))

  do i = 1, nedgc

    e = pedge(1,i)
    ee = pedge(2,i)

    if ( fvl(facn,e) == fl ) then
      a = e
    else if ( fvl(facn,ee) == fl ) then
      a = ee
    else
      cycle
    end if

    ca = pedge(3,i) / 10
    cb = mod(pedge(3,i),10)

    if ( ca == 2 ) then

      la = fvl(loc,a)

      if ( cb == 2 ) then

        lb = fvl(loc,fvl(succ,a))
        s = (vcl(k,la) - vcl(k,lw)) / dir(k)
        t = (vcl(k,lb) - vcl(k,lw)) / dir(k)

        if ( 0.0D+00 < s ) then

          if ( min ( s, t ) < tmin .or. imin == 0 ) then

            if ( s < t ) then

              if ( ccwfl < 0 ) then
                imin = a
                lw1 = la
                tmin = s
              end if

            else

              if ( 0 < ccwfl ) then
                imin = a
                lw1 = lb
                tmin = t
              end if

            end if

          end if

        end if

      else

        l = fvl(loc,e)

        if ( ( l == la .and. ccwfl < 0 ) .or. &
             ( l /= la .and. 0 < ccwfl ) ) then

          t = ( vcl(k,l) - vcl(k,lw) ) / dir(k)

          if ( ntol < t ) then

            if ( t < tmin .or. imin == 0 ) then
              lw1 = l
              imin = a
              tmin = t
            end if

          end if

        end if

      end if

    else if ( cb == 2 ) then

      la = fvl(loc,a)
      l = fvl(loc,fvl(succ,e))

      if ( ( l == la .and. ccwfl < 0 ) .or. &
           ( l /= la .and. 0 < ccwfl ) ) then

        t = ( vcl(k,l) - vcl(k,lw) ) / dir(k)

        if ( ntol < t ) then

          if ( t < tmin .or. imin == 0 ) then
            lw1 = l
            imin = a
            tmin = t
          end if

        end if

      end if

    else

      la = fvl(loc,e)
      lb = fvl(loc,fvl(succ,e))

      dir1(1:3) = vcl(1:3,la) - vcl(1:3,lb)
      rhs(1:3) = vcl(1:3,la) - vcl(1:3,lw)

      cp(1) = dir(2)*dir1(3) - dir(3)*dir1(2)
      cp(2) = dir(3)*dir1(1) - dir(1)*dir1(3)
      cp(3) = dir(1)*dir1(2) - dir(2)*dir1(1)

      l = 1
      if ( abs(cp(1)) < abs(cp(2)) ) then
        l = 2
      end if

      if ( abs(cp(l)) < abs(cp(3)) ) then
        l = 3
      end if

      if ( l == 1 ) then
        t = (rhs(2)*dir1(3) - rhs(3)*dir1(2))/cp(1)
      else if ( l == 2 ) then
        t = (rhs(3)*dir1(1) - rhs(1)*dir1(3))/cp(2)
      else
        t = ( rhs(1) * dir1(2) - rhs(2) * dir1(1) ) / cp(3)
      end if

      if ( ntol < t ) then

        if ( t < tmin .or. imin == 0 ) then
          imin = -a
          tmin = t
          isave = i
        end if

      end if

    end if

  end do

  if ( imin == 0 ) then
    return
  end if

  if ( imin < 0 ) then

    if ( nvc < lv ) then

      if ( -imin == cedge(2,0) ) then
        lw1 = lv
        go to 90
      end if

    end if

    n = n + 1

    if ( maxvc < n ) then
      ierr = 14
      return
    end if

    lw1 = n
    vcl(1:3,n) = vcl(1:3,lw) + tmin * dir(1:3)

90  continue

    if ( abs ( facep(2,fl) ) /= abs ( facep(3,fl) ) ) then
      do i = 1, 3
        j = pedge(i,isave)
        pedge(i,isave) = pedge(i,nedgc)
        pedge(i,nedgc) = j
      end do
      nedgc = nedgc - 1
    end if

  end if

  nce = nce + 1
  cedge(1,nce) = lw1
  cedge(2,nce) = -abs(imin)
!
!  If vertex of cut polygon has appeared before, then cut polygon
!  is simply-connected (non-simple), so reject cut plane.
!
110 continue

  if ( lw1 == lv ) then
    go to 150
  end if

  if ( lw1 <= nvc ) then

    do i = 1, nev
      if ( lw1 == ev(i) ) then
        if ( msglvl == 4 ) then
          write ( *,600) 'rejected due to simply-connected polygon: case 1'
        end if
        return
      end if
    end do

    if ( maxev <= nev ) then
      ierr = 328
      return
    end if

    nev = nev + 1
    ev(nev) = lw1

  else

    do i = nvc+1, n-1

      do j = 1, 3
        cmax = max ( abs(vcl(j,i)), abs(vcl(j,n)) )
        if ( tol * cmax < abs(vcl(j,i) - vcl(j,n)) .and. &
          tol < cmax ) then
          go to 140
        end if
      end do

      if ( msglvl == 4 ) then
        write ( *,600) 'rejected due to simply-connected polygon: case 2'
      end if

      return

140   continue

    end do

  end if
!
!  Compute dihedral angles due to cut plane at edge. If any angle
!  is too small, reject cut plane.
!
150 continue

  if ( fr == 0 ) then

    f = fl
    dof = ( abs(facep(2,f)) == abs(facep(3,f)) )
    eflag = ( abs(facep(2,f) ) /= p .or. ( dof .and. ccwfl < 0 ) )

  else

    f = fr
    dof = ( abs(facep(2,f)) == abs(facep(3,f)) )
    sf = ccwfl

    if ( dof ) then

      e = cedge(2,nce)
      la = fvl(loc,e)
      lb = fvl(loc,fvl(succ,e))

      if ( 0 < (lb - la) * ccwfl ) then
        e = fvl(edgc,e)
      else
        e = fvl(edga,e)
      end if

      l = fvl(loc,e)
      if ( l == la ) then
        sf = -ccwfl
      end if

    end if

    eflag = ( abs(facep(2,f)) /= p .or. ( dof .and. sf < 0 ) )

  end if

  dotp = - dot_product ( nrmlc(1:3), nrml(1:3,f) )

  if ( 1.0D+00 - tol < abs ( dotp ) ) then
    dotp = sign ( 1.0D+00, dotp )
  end if

  if ( eflag ) then
    dotp = -dotp
  end if

  angr = pi - acos(dotp)
  dir1(1) = nrmlc(2)*dir(3) - nrmlc(3)*dir(2)
  dir1(2) = nrmlc(3)*dir(1) - nrmlc(1)*dir(3)
  dir1(3) = nrmlc(1)*dir(2) - nrmlc(2)*dir(1)

  dotp = dot_product ( dir1(1:3), nrml(1:3,f) )

  if ( eflag ) then
    dotp = -dotp
  end if

  if ( 0.0D+00 < dotp ) then
    angr = pi2 - angr
  end if

  if ( fr == 0 ) then

    ang = pi

  else

    a = cedge(2,nce)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))

    if ( 0 < (lb - la) * ccwfl ) then
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

  if ( lw1 /= lv ) then
    go to 10
  end if
!
!  Determine if cut polygon is outer or inner boundary by summing
!  the exterior angles (which lie in range (-PI,PI)). A sum of 2*PI
!  (-2*PI) means that polygon is outer (inner). Cut polygon is
!  rejected in the latter case.
!
  s = 0.0D+00
  la = cedge(1,nce-1)
  lb = cedge(1,0)
  dir1(1) = vcl(1,lb) - vcl(1,la)
  dir1(2) = vcl(2,lb) - vcl(2,la)
  dir1(3) = vcl(3,lb) - vcl(3,la)
  dir1sq = dir1(1)**2 + dir1(2)**2 + dir1(3)**2

  do i = 0, nce-1

    dir(1:3) = dir1(1:3)

    dirsq = dir1sq
    la = lb
    lb = cedge(1,i+1)

    dir1(1:3) = vcl(1:3,lb) - vcl(1:3,la)

    dir1sq = dir1(1)**2 + dir1(2)**2 + dir1(3)**2
    dotp = (dir(1)*dir1(1) + dir(2)*dir1(2) + dir(3)*dir1(3))/ &
      sqrt(dirsq*dir1sq)

    if ( 1.0D+00 - tol < abs ( dotp ) ) then
      dotp = sign ( 1.0D+00, dotp )
    end if

    ang = acos(dotp)

    if ( kmax == 1 ) then
      cp(1) = dir(2)*dir1(3) - dir(3)*dir1(2)
    else if ( kmax == 2 ) then
      cp(2) = dir(3)*dir1(1) - dir(1)*dir1(3)
    else
      cp(3) = dir(1)*dir1(2) - dir(2)*dir1(1)
    end if

    if ( cp(kmax)*nrmlc(kmax) < 0.0D+00 ) then
      ang = -ang
    end if

    s = s + ang

  end do

  if ( s < 0.0D+00 ) then
    if ( msglvl == 4 ) then
      write ( *, '(a)' ) 'Rejected due to inner boundary'
      return
    end if
  end if
!
!  Move edges incident on LV (if <= NVC), EV(1:NEV) to end of PEDGE.
!
  if ( lv <= nvc ) then
    l = lv
    ee = 0
  else
    ee = 1
  end if

  do e = ee, nev

    if ( 0 < e ) then
      l = ev(e)
    end if

    k = 1

    do while ( k <= nedgc )

      a = pedge(1,k)
      la = fvl(loc,a)
      lb = fvl(loc,fvl(succ,a))

      if ( la == l .or. lb == l ) then
        do i = 1, 3
          j = pedge(i,k)
          pedge(i,k) = pedge(i,nedgc)
          pedge(i,nedgc) = j
        end do
        nedgc = nedgc - 1
      else
        k = k + 1
      end if

    end do

  end do
!
!  Determine if cut face contains any inner polygons by checking
!  if the remaining edges of PEDGE intersect interior of cut face.
!
  do i = 1, nedgc

    a = pedge(1,i)
    ca = pedge(3,i)/10
    cb = mod(pedge(3,i),10)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))

    if ( ca == 2 ) then
      cp(1) = vcl(1,la)
      cp(2) = vcl(2,la)
      cp(3) = vcl(3,la)
    else if ( cb == 2 ) then
      cp(1) = vcl(1,lb)
      cp(2) = vcl(2,lb)
      cp(3) = vcl(3,lb)
    else
      dir(1) = vcl(1,lb) - vcl(1,la)
      dir(2) = vcl(2,lb) - vcl(2,la)
      dir(3) = vcl(3,lb) - vcl(3,la)
      t = (nrmlc(4) - nrmlc(1)*vcl(1,la) - nrmlc(2)*vcl(2,la) - &
           nrmlc(3)*vcl(3,la))/(nrmlc(1)*dir(1) + nrmlc(2)*dir(2) &
           + nrmlc(3)*dir(3))
      cp(1) = vcl(1,la) + t*dir(1)
      cp(2) = vcl(2,la) + t*dir(2)
      cp(3) = vcl(3,la) + t*dir(3)
    end if

    call ptpolg(3,3,nce,2,cedge,vcl,cp,nrmlc,dtol,inout)

    if ( inout == 1 ) then
      if ( msglvl == 4 ) then
        write ( *,600) 'rejected due to hole polygon'
      end if
      return
    end if

  end do

  rflag = .true.

  if ( msglvl == 4 ) then
    write ( *,'(a)' ) 'cedge(1:2), cdang'
    do i = 1, nce
      write ( *,610) i,cedge(1,i),cedge(2,i),cdang(i)*180.0D+00 / pi
    end do
  end if

  600 format (4x,a)
  610 format (4x,3i7,f12.5)

  return
end
