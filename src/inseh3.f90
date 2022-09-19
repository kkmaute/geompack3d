subroutine inseh3 ( a, b, nvert, maxfv, facep, fvl, eang, ierr )

!*****************************************************************************80
!
!! INSEH3 inserts an edge into the polyhedral decomposition data structure.
!
!  Discussion:
!
!    This routine inserts an edge on a face of polyhedral decomposition data
!    structure that joins a hole on face to outer boundary polygon
!    of face.
!
!    It is assumed that the edge is entirely inside the face.
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
!    Input, integer ( kind = 4 ) A, B, the indices of FVL for vertices on same face, one
!    on hole and the other on outer boundary.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfv

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) nvert
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
  fvl(facn,i) = f
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
  sp = facep(2,f)
  sq = facep(3,f)

  if ( sq == 0 ) then

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

  return
end
