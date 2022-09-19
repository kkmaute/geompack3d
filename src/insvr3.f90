subroutine insvr3 ( a, nvc, nvert, maxfv, vcl, fvl, eang, ierr )

!*****************************************************************************80
!
!! INSVR3 inserts a point into the polyhedral decomposition database.
!
!  Discussion:
!
!    This routine inserts a vertex on an edge of polyhedral decomposition
!    data structure.
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
!    Input, integer ( kind = 4 ) A, the index of FVL specifying edge containing
!    inserted vertex.
!
!    Input/output, integer ( kind = 4 ) NVC, the number of vertex coordinates.
!
!    Input/output, integer ( kind = 4 ) NVERT, the number of positions used in FVL,
!    EANG arrays.
!
!    Input, integer ( kind = 4 ) MAXFV, the maximum size available for FVL, EANG arrays.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:NVC+1), the vertex coordinate list;
!    VCL(*,NVC+1) are coordinates of vertex to be inserted.
!
!    Input/output, integer ( kind = 4 ) FVL(1:6,1:NVERT), the face vertex list.
!
!    Input/output, real ( kind = 8 ) EANG(1:NVERT), the edge angles.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxfv
  integer ( kind = 4 ) nvc

  integer ( kind = 4 ) a
  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angnxt
  integer ( kind = 4 ) b
  logical              bflag
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  real    ( kind = 8 ) eang(maxfv)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(6,maxfv)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) la
  integer ( kind = 4 ) lb
  integer ( kind = 4 ) li
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lnxt
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) vcl(3,nvc)

  ierr = 0
  nvc = nvc + 1
  b = fvl(succ,a)
  la = fvl(loc,a)
  lb = fvl(loc,b)
  i = a
!
!  Find start edge of FVL if AB lies on boundary of decomposition.
!
  c = i

  do

    d = fvl(edga,c)

    if ( d == 0 ) then
      exit
    end if

    if ( d == i ) then
      exit
    end if

    c = d

  end do

  bflag = ( d == 0 )
!
!  Insert new entry of FVL for each face containing edge AB.
!
  i = c
  k = nvert

  do

    j = fvl(succ,i)
    k = k + 1

    if ( maxfv < k ) then
      ierr = 15
      return
    end if

    fvl(loc,k) = nvc
    fvl(facn,k) = fvl(facn,i)
    fvl(succ,k) = j
    fvl(pred,k) = i
    fvl(edga,k) = 0
    fvl(edgc,k) = 0
    fvl(succ,i) = k
    fvl(pred,j) = k
    eang(k) = -1.0D+00
    i = fvl(edgc,i)

    if ( i == 0 .or. i == c ) then
      exit
    end if

  end do

  nvert = k
!
!  Set EDGA, EDGC, and EANG fields.
!
  i = c
  l = fvl(edgc,i)
  ang = eang(i)

30 continue

  k = fvl(succ,i)
  j = fvl(succ,k)
  li = fvl(loc,i)
  lj = fvl(loc,j)
  n = fvl(succ,l)
  lnxt = fvl(edgc,l)
  angnxt = eang(l)

  if ( li < lj ) then

    if ( fvl(loc,l) == li ) then
      fvl(edga,k) = n
      fvl(edgc,n) = k
      eang(n) = ang
    else
      fvl(edgc,i) = n
      fvl(edga,n) = i
      fvl(edga,k) = l
      fvl(edgc,l) = k
      eang(l) = ang
      if ( lnxt == 0 ) then
        fvl(edga,l) = 0
      end if
    end if

  else

    if ( fvl(loc,l) == li ) then
      fvl(edgc,k) = n
      fvl(edga,n) = k
      eang(k) = ang
      fvl(edga,i) = l
      fvl(edgc,l) = i
      eang(l) = ang
      if ( lnxt == 0 ) then
        fvl(edga,l) = 0
      end if
    else
      fvl(edgc,k) = l
      fvl(edga,l) = k
      eang(k) = ang
      fvl(edga,i) = n
      fvl(edgc,n) = i
      eang(n) = ang
    end if

    if ( bflag .and. i == c ) then
      fvl(edgc,i) = 0
      eang(i) = -1.0D+00
    end if

  end if

  i = l
  l = lnxt
  ang = angnxt
  if ( l /= 0 .and. i /= c ) go to 30

  return
end
