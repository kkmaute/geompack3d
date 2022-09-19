subroutine sfc1mf ( p, cntr, mean, nf, indf, meanf, angacc, mxcos, dtol, nvc, &
  maxvc, vcl, facep, nrml, fvl, eang, nrmlc, nce, cedge, cdang, aflag, &
  iwk, ierr )

!*****************************************************************************80
!
!! SFC1MF seeks a separator or cut face in a convex polyhedron.
!
!  Discussion:
!
!    This routine attempts to find a separator or cut face in a convex
!    polyhedron P based on a mesh distribution function by starting with
!    an edge of a polyhedron with a large enough dihedral angle and big
!    difference between mean mdf in 2 faces incident on edge.
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
!    Input, real ( kind = 8 ) MEAN, the mean MDF value in polyhedron.
!
!    Input, integer ( kind = 4 ) NF, the number of faces in polyhedron.
!
!    Input, integer ( kind = 4 ) INDF(1:NF), the indices in FACEP of faces of polyhedron.
!
!    Input, real ( kind = 8 ) MEANF(1:NF), the mean mdf value associated
!    with faces of polyhedron.
!
!    Input, real ( kind = 8 ) ANGACC, the min acceptable dihedral angle
!    in radians produced by a cut face.
!
!    Input, real ( kind = 8 ) MXCOS, the max cosine of angle allowed for
!    angles subtended by new subedges with respect to centroid.
!
!    Input, real ( kind = 8 ) DTOL, the absolute tolerance to determine if
!    point is on cut plane.
!
!    Input, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input, integer ( kind = 4 ) MAXVC, the maximum size available for VCL array.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate
!    list.  On output, some temporary or permanent entries may be added at
!    end of array.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, real ( kind = 8 ) NRML(1:3,1:*), the unit normal vectors for faces.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, real ( kind = 8 ) EANG(1:*), the edge angles.
!
!    Output, real ( kind = 8 ) NRMLC(1:4), the unit normal vector of cut
!    plane plus right hand side constant term of plane equation (if
!    acceptable).
!
!    Output, integer ( kind = 4 ) NCE, the number of edges in cut face (if acceptable); it is
!    assumed there is enough space in the following two arrays.
!
!    Output, integer ( kind = 4 ) CEDGE(1:2,0:NCE), real ( kind = 8 ) CDANG(1:NCE), the
!    information describing cut polygon as output by routine SEPFAC.
!
!    Output, logical AFLAG, is TRUE iff separator face is found and acceptable.
!
!    Workspace, integer IWK(1:NF).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxvc
  integer ( kind = 4 ) nf

  integer ( kind = 4 ) a
  logical              aflag
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angacc
  integer ( kind = 4 ) b
  real    ( kind = 8 ) cdang(*)
  real    ( kind = 8 ) ce
  integer ( kind = 4 ) cedge(2,0:*)
  real    ( kind = 8 ) cn
  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) diffi
  real    ( kind = 8 ) dtol
  integer ( kind = 4 ) e(2)
  real    ( kind = 8 ) eang(*)
  real    ( kind = 8 ) edg(3)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) enr(3)
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  real    ( kind = 8 ), parameter, dimension ( 3 ) :: fract = (/ &
    0.5D+00, 0.4D+00, 0.6D+00 /)
  integer ( kind = 4 ) fvl(6,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) indf(nf)
  integer ( kind = 4 ) iwk(nf)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  real    ( kind = 8 ) leng
  integer ( kind = 4 ), parameter :: loc = 1
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) meanf(nf)
  real    ( kind = 8 ) mndang
  integer ( kind = 4 ), parameter :: msglvl = 0
  real    ( kind = 8 ) mxcos
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nce
  logical              neg
  real    ( kind = 8 ) nrml(3,*)
  real    ( kind = 8 ) nrmlc(4)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) p
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) r(2)
  real    ( kind = 8 ) ratio
  integer ( kind = 4 ) sp
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(3,maxvc)
!
!  Find up to 2 candidates for starting edge of cut face.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  aflag = .false.
  mndang = min ( pi * 0.9D+00, angacc * 3.0D+00 )
  n = 0

  do i = 1, nf
    f = indf(i)
    iwk(i) = facep(1,f)
    facep(1,f) = i
  end do

  do i = 1, nf

    diffi = meanf(i) - mean
    f = indf(i)

    if ( abs ( facep(2,f) ) == p ) then
      sp = facep(2,f)
    else
      sp = facep(3,f)
    end if

    a = iwk(i)
    la = fvl(loc,a)

20  continue

    b = fvl(succ,a)
    lb = fvl(loc,b)

    if ( 0 < ( lb - la ) * sp ) then

      if ( eang(a) < mndang ) then
        go to 30
      end if

      f = fvl(facn,fvl(edgc,a))
      j = facep(1,f)

      if ( 0.0D+00 <= diffi * ( meanf(j) - mean ) ) then
        go to 30
      end if

      if ( 0.0D+00 < diffi ) then
        ratio = meanf(i) / meanf(j)
      else
        ratio = meanf(j) / meanf(i)
      end if

      if ( ratio < 1.5D+00 ) then
        go to 30
      end if

      if ( n <= 1 ) then

        n = n + 1
        e(n) = a
        r(n) = ratio

      else

        if ( r(1) < r(2) ) then
          k = 1
        else
          k = 2
        end if

        if ( r(k) < ratio ) then
          e(k) = a
          r(k) = ratio
        end if

      end if

    end if

30  continue

    a = b
    la = lb

    if ( a /= iwk(i) ) then
      go to 20
    end if

  end do

  do i = 1, nf
    facep(1,indf(i)) = iwk(i)
  end do

  if ( n == 2 ) then
    if ( r(1) < r(2) ) then
      j = e(1)
      e(1) = e(2)
      e(2) = j
    end if
  end if
!
!  For each candidate edge, try 3 different cut planes.
!
  do k = 1, n

    a = e(k)
    la = fvl(loc,a)
    lb = fvl(loc,fvl(succ,a))

    if ( lb < la ) then
      j = la
      la = lb
      lb = j
    end if

    cedge(1,0) = lb
    cedge(1,1) = la
    cedge(2,1) = a

    edg(1:3) = vcl(1:3,lb) - vcl(1:3,la)

    leng = sqrt ( edg(1)**2 + edg(2)**2 + edg(3)**2 )
    edg(1:3) = edg(1:3) / leng
    f = fvl(facn,fvl(edgc,a))
    enr(1) = edg(2) * nrml(3,f) - edg(3) * nrml(2,f)
    enr(2) = edg(3) * nrml(1,f) - edg(1) * nrml(3,f)
    enr(3) = edg(1) * nrml(2,f) - edg(2) * nrml(1,f)
    neg = ( abs ( facep(2,f) ) /= p )

    if ( neg ) then
      enr(1:3) = -enr(1:3)
    end if

    do i = 1, 3

      ang = eang(a) * fract(i)

      if ( msglvl == 4 ) then
        write ( *,600) a,lb,la,f,p, eang(a)*180.0D+00/ pi, &
          ang * 180.0D+00 / pi
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

  600 format (' sfc1mf: a,lb,la,f,p,eang(a),ang =',5i5,2f9.4)

  return
end
