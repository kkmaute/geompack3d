subroutine insed3 ( a, b, nface, nvert, npf, maxfp, maxfv, maxpf, facep, &
  factyp, nrml, fvl, eang, hfl, pfl, ierr )

!*****************************************************************************80
!
!! INSED3 inserts an edge into the polyhedral decomposition data structure.
!
!  Discussion:
!
!    This routine inserts an edge on a face of polyhedral decomposition data
!    structure.  It is assumed that the edge is entirely inside the face.
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
!    Input/output, integer ( kind = 4 ) A, B, indices of FVL for nonadjacent vertices
!    on same face.  On output, A is the index in FVL of the new edge; LOC
!    field of output A is the same as that of the input A.
!
!    Input/output, integer ( kind = 4 ) NFACE, the number of faces or positions used
!    in FACEP array.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input/output, integer ( kind = 4 ) NPF, the number of positions used in PFL array.
!
!    Input, integer ( kind = 4 ) MAXFP, the maximum size available for FACEP, FACTYP,
!    NRML arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXPF, the maximum size available for PFL array.
!
!    Input/output, integer ( kind = 4 ) FACEP(1:3,1:NFACE), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FACTYP(1:NFACE), the face types.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), the unit normal
!    vectors for faces; outward normal corresponds to counterclockwise
!    traversal of face from polyhedron with index FACEP(2,F).
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Input, integer ( kind = 4 ) HFL(1:*), the head pointer to face indices in PFL
!    for each polyhedron.
!
!    Input/output, integer ( kind = 4 ) PFL(1:2,1:NPF), the list of signed face
!    indices for each polyhedron.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfp
  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) maxpf

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,maxfp)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) factyp(maxfp)
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hfl(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nface
  real    ( kind = 8 ) nrml(3,maxfp)
  integer ( kind = 4 ) npf
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) pfl(2,maxpf)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) sp
  integer ( kind = 4 ) sq
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  f = fvl(facn,a)
  i = nvert + 1
  j = nvert + 2
  nvert = j

  if ( maxfv < nvert ) then
    ierr = 15
    return
  end if

  fvl(loc,i) = fvl(loc,a)
  fvl(succ,i) = b
  k = fvl(pred,a)
  fvl(pred,i) = k
  fvl(succ,k) = i
  fvl(edga,i) = j
  fvl(edgc,i) = j
  fvl(loc,j) = fvl(loc,b)
  fvl(facn,j) = f
  fvl(succ,j) = a
  k = fvl(pred,b)
  fvl(pred,j) = k
  fvl(succ,k) = j
  fvl(edga,j) = i
  fvl(edgc,j) = i
  fvl(pred,a) = j
  fvl(pred,b) = i
  eang(i) = pi
  eang(j) = pi
  nface = nface + 1

  if ( maxfp < nface ) then
    ierr = 16
    return
  end if

  facep(1,f) = a
  sp = facep(2,f)
  sq = facep(3,f)
  facep(1,nface) = b
  facep(2,nface) = sp
  facep(3,nface) = sq
  factyp(nface) = factyp(f)
  nrml(1,nface) = nrml(1,f)
  nrml(2,nface) = nrml(2,f)
  nrml(3,nface) = nrml(3,f)
  k = b

  do

    fvl(facn,k) = nface
    k = fvl(succ,k)

    if ( k == b ) then
      exit
    end if

  end do

  g = hfl( abs ( sp ) )
  npf = npf + 1

  if ( maxpf < npf ) then
    ierr = 17
    return
  end if

  pfl(1,npf) = sign ( nface, sp )
  pfl(2,npf) = pfl(2,g)
  pfl(2,g) = npf

  if ( sq /= 0 ) then
    g = hfl( abs ( sq ) )
    npf = npf + 1

    if ( maxpf < npf ) then
      ierr = 17
      return
    end if

    pfl(1,npf) = sign ( nface, sq )
    pfl(2,npf) = pfl(2,g)
    pfl(2,g) = npf
  else
    if ( 0 < (fvl(loc,b) - fvl(loc,a))*sp ) then
      fvl(edga,i) = 0
      fvl(edgc,j) = 0
      eang(j) = -1.0D+00
    else
      fvl(edga,j) = 0
      fvl(edgc,i) = 0
      eang(i) = -1.0D+00
    end if
  end if

  a = i

  return
end
