subroutine xpghpl ( ar, br, dr, g, maxsv, k, head, svcl, sfvl, empty, ierr )

!*****************************************************************************80
!
!! XPGHPL intersects a convex polygon with a half plane.
!
!  Discussion:
!
!    This routine determines the intersection of a convex polygon in SVCL, SFVL
!    data structure (of routine SHRNK3) with half-plane determined
!    by line G.
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
!    Input, real ( kind = 8 ) AR, BR, DR, the half-plane equation is
!    AR*X + BR*Y <= DR.
!
!    Input, integer ( kind = 4 ) G, the index of line (actually projection of plane).
!
!    Input, integer ( kind = 4 ) MAXSV, the maximum size available for SVCL, SFVL arrays.
!
!    Input/output, integer ( kind = 4 ) K, the index of last entry used in SVCL, SFVL arrays.
!
!    Input/output, integer ( kind = 4 ) HEAD, the head pointer (index of SFVL) to vertex of
!    convex polygon.
!
!    Input/output, real ( kind = 8 ) SVCL(1:3,1:MAXSV),
!    integer ( kind = 4 ) SFVL(1:5,1:MAXSV); see routine SHRNK3.
!
!    Output, logical EMPTY, is TRUE if intersection is empty or degenerate;
!    else FALSE.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxsv

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ar
  integer ( kind = 4 ) b
  real    ( kind = 8 ) br
  real    ( kind = 8 ) dn
  real    ( kind = 8 ) dp
  real    ( kind = 8 ) dr
  real    ( kind = 8 ) dv
  real    ( kind = 8 ) dx
  real    ( kind = 8 ) dy
  integer ( kind = 4 ), parameter :: edgv = 5
  logical              empty
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fdel
  integer ( kind = 4 ) fkeep
  integer ( kind = 4 ) g
  integer ( kind = 4 ) head
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) ionl
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jj
  integer ( kind = 4 ) k
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) sfvl(5,maxsv)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) svcl(3,maxsv)
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  empty = .false.
  f = sfvl(facn,head)
  fdel = 0
  fkeep = 0
  ionl = 0
  dv = abs(dr)*tol
  dn = dr - dv
  dp = dr + dv

  if ( abs ( dr ) <= tol ) then
    dp = tol + tol
    dn = -dp
  end if
!
!  Determine if intersection is empty, entire polygon, or neither.
!
  i = head

10 continue

  j = sfvl(loc,i)
  dv = ar * svcl(1,j) + br * svcl(2,j)

  if ( dv < dn ) then

    if ( fkeep == 0 ) then
      fkeep = i
    end if

  else if ( dp < dv ) then

    if ( fdel == 0 ) then
      fdel = i
    end if

    sfvl(loc,i) = -j
    sfvl(facn,i) = 0
!
!  A deleted vertex has a negative LOC value and 0 FACN value.
!
   else
      ionl = i
      sfvl(facn,i) = -f
!
!  A vertex lying on line has a negative FACN value.
!
   end if

  i = sfvl(succ,i)
  if ( i /= head ) go to 10

  if ( fdel == 0 ) then
    if( ionl /= 0 ) then
      sfvl(facn,ionl) = f
      i = sfvl(succ,ionl)
      sfvl(facn,i) = f
      i = sfvl(pred,ionl)
      sfvl(facn,i) = f
    end if
    return
  else if ( fkeep == 0 ) then
    empty = .true.
    head = 0
    return
  end if
!
!  Update polygon to include part of line and remove at least 1 vertex.
!
  i = fdel

20 continue

  i = sfvl(pred,i)
  j = sfvl(loc,i)
  if ( j < 0 ) go to 20

  if ( sfvl(facn,i) < 0 ) then

    a = i
    sfvl(facn,i) = f

  else

    ii = sfvl(succ,i)
    jj = -sfvl(loc,ii)
    dx = svcl(1,jj) - svcl(1,j)
    dy = svcl(2,jj) - svcl(2,j)
    t = ( dr - ar * svcl(1,j) - br * svcl(2,j) ) / ( ar * dx + br * dy )

    if ( t <= tol .or. 1.0D+00 - tol <= t ) then
      ierr = 313
      return
    end if

    if ( ii /= fdel ) then

      a = ii

    else

      k = k + 1

      if ( maxsv < k ) then
        ierr = 13
        return
      end if

      a = k
      sfvl(pred,k) = i
      sfvl(succ,i) = k
      jj = k

    end if

    svcl(1,jj) = svcl(1,j) + t * dx
    svcl(2,jj) = svcl(2,j) + t * dy
    svcl(3,jj) = 0.0D+00
    sfvl(loc,a) = jj
    sfvl(facn,a) = f

  end if

  i = fdel

  do

    i = sfvl(succ,i)
    j = sfvl(loc,i)

    if ( 0 <= j ) then
      exit
    end if

  end do

  if ( sfvl(facn,i) < 0 ) then

    b = i
    sfvl(facn,i) = f

  else

    ii = sfvl(pred,i)
    jj = -sfvl(loc,ii)
    b = ii
    dx = svcl(1,jj) - svcl(1,j)
    dy = svcl(2,jj) - svcl(2,j)
    t = ( dr - ar * svcl(1,j) - br * svcl(2,j) ) / ( ar * dx + br * dy )

    if ( t <= tol .or. 1.0D+00 - tol <= t ) then
      ierr = 313
      return
    end if

    svcl(1,jj) = svcl(1,j) + t * dx
    svcl(2,jj) = svcl(2,j) + t * dy
    sfvl(loc,b) = jj
    sfvl(facn,b) = f

  end if

  sfvl(succ,a) = b
  sfvl(pred,b) = a
  sfvl(edgv,a) = -g
  head = a

  return
end
