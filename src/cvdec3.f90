subroutine cvdec3 ( angacc, rdacc, nvc, nface, nvert, npolh, npf, maxvc, &
  maxfp, maxfv, maxhf, maxpf, maxiw, maxwk, vcl, facep, factyp, nrml, fvl, &
  eang, hfl, pfl, iwk, wk, ierr )

!*****************************************************************************80
!
!! CVDEC3 decomposes polyhedra into convex parts.
!
!  Discussion:
!
!    This routine is given one or more polyhedra in polyhedral decomposition
!    data structure, and decomposes the polyhedra into convex parts.
!
!    It is assumed all faces are simple (any faces with holes should be
!    decomposed into simple polygons), each face appears at most
!    once in a polyhedron (double-occurring faces are not allowed),
!    and no interior holes occur in any polyhedra.
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
!    Input, real ( kind = 8 ) ANGACC, the minimum acceptable dihedral angle
!    in radians produced by cut faces.
!
!    Input, real ( kind = 8 ) RDACC, the minimum acceptable relative
!    distance between cut planes and vertices not on plane.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates or
!    positions used in VCL.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used
!    in FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPOLH, the number of polyhedra or positions
!    used in HFL array.
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
!    Input, integer ( kind = 4 ) MAXIW, the maximum size available for IWK array;
!    should be about 5 times number of edges in any polyhedron.
!
!    Input, integer ( kind = 4 ) MAXWK, the maximum size available for WK array; should
!    be about number of edges in any polyhedron.
!
!    Input/output, real ( kind = 8 ) VCL(1:3,1:NVC), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list: row 1 is
!    head pointer, rows 2 and 3 are signed polyhedron indices.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types: useful for
!    specifying types of boundary faces; entries must be greater than or
!    equal to 0; any new interior aces (not part of previous face) has
!    face type set to 0.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces; outward normal corresponds to counterclockwise
!    traversal of face from polyhedron with index |FACEP(2,F)|.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list; see
!    routine DSPHDC.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the angles at edges
!    common to 2 faces in a polyhedron; EANG(J) corresponds to FVL(*,J),
!    determined by EDGC field.
!
!    Input/output, integer ( kind = 4 ) HFL(1:NPOLH), the head pointer to face indices
!    in PFL for each polyhedron.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:NPF), list of signed face indices for
!    each polyhedron; row 2 used for link.
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

  real    ( kind = 8 ) angacc
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) hfl(maxhf)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iwk(maxiw)
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nface
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) npolh
  real    ( kind = 8 ) nrml(3,maxfp)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) piptol
  integer ( kind = 4 ) pfl(2,maxpf)
  integer ( kind = 4 ), parameter :: pred = 4
  real    ( kind = 8 ) rdacc
  logical              rflag
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) u
  real    ( kind = 8 ) vcl(3,maxvc)
  real    ( kind = 8 ) wk(maxwk)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  piptol = pi + tol

10 continue

  l = 0
  n = 0
  u = 1

  do

    if ( piptol < eang(u) ) then

      call resedg(u,angacc,rdacc,nvc,nface,nvert,npolh,npf,maxvc, &
        maxfp,maxfv,maxhf,maxpf,maxiw/3,maxwk,vcl,facep,factyp, &
        nrml,fvl,eang,hfl,pfl,rflag,iwk,wk, ierr )

      if ( ierr /= 0 ) then
        return
      end if

      if ( rflag ) then
        n = n + 1
      else
        if ( l == 0 ) then
          l = u
        end if
      end if

    end if

    u = u + 1

    if ( nvert < u ) then
      exit
    end if

  end do

  if ( 0 < l ) then
    if ( n == 0 ) then
      ierr = 327
      return
    else
      go to 10
    end if
  else
    if ( 0 < n ) then
      go to 10
    end if
  end if

  return
end
