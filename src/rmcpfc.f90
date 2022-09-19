subroutine rmcpfc ( nface, nvert, hvl, nrml, fvl, eang, iwk )

!*****************************************************************************80
!
!! RMCPFC removes coplanar adjacent polyhedron faces from the data base.
!
!  Discussion:
!
!    This routine removes coplanar adjacent faces from a convex polyhedron
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
!    Input/output, integer ( kind = 4 ) NFACE, number of faces in convex polyhedron.
!
!    Input/output, integer ( kind = 4 ) NVERT, size of FVL, EANG arrays; must be
!    greater than or equal to twice the number of edges of polyhedron.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NFACE), head vertex list.  On output,
!    first NFACE entries indicate resulting faces.
!
!    Input/output, real ( kind = 8 ) NRML(1:3,1:NFACE), unit outward normals
!    of faces.
!
!    Input/output, integer ( kind = 4 ) FVL(1:5,1:NVERT), face vertex list; see routine
!    DSCPH.  On output, may contain some unused columns, indicated by
!    nonpositive LOC values.
!
!    Input, real ( kind = 8 ) EANG(1:NVERT), angles at edges common to 2
!    faces; EANG(I) corresponds to FVL(*,I).
!
!    Workspace, integer IWK(1:NVERT).
!
  implicit none

  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  real    ( kind = 8 ) eang(nvert)
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,nvert)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) iwk(nvert)
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) np
  real    ( kind = 8 ) nrml(3,nface)
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  real    ( kind = 8 ) pitol
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) tol

  tol = 100.0D+00 * epsilon ( tol )
  pitol = pi - tol
  np = 0

  do i = 1, nvert
    iwk(i) = 0
    if ( fvl(loc,i) <= 0) iwk(i) = 1
  end do

  do ii = 1, nvert

    if ( iwk(ii) == 1 ) then
      cycle
    end if

    i = ii
    a = 0
    b = 0

20  continue

    iwk(i) = 1

    if ( pitol <= eang(i) ) then
      np = np + 1
      fvl(loc,i) = 0
      a = i
    else
      b = i
    end if

    i = fvl(succ,fvl(edgv,i))
    if ( i /= ii) go to 20

    if ( a == 0 .or. b == 0 ) then
      cycle
    end if

    a = b

30  continue

    i = fvl(edgv,a)
    j = fvl(succ,i)

    if ( pitol <= eang(j) ) then

40    continue

      j = fvl(succ,fvl(edgv,j))

      if ( pitol <= eang(j) ) then
        go to 40
      end if

      fvl(succ,i) = j
      fvl(pred,j) = i

    end if

    a = j
    if ( a /= b) go to 30

  end do

  if ( np == 0 ) then
    return
  end if

  do i = 1, nvert

    if ( iwk(i) == 2 .or. fvl(loc,i) <= 0)then
      cycle
    end if

    f = fvl(facn,i)
    hvl(f) = i
    j = i

    do

      iwk(j) = 2
      g = fvl(facn,j)

      if ( g /= f ) then
        fvl(facn,j) = f
        hvl(g) = 0
      end if

      j = fvl(succ,j)

      if ( j == i ) then
        exit
      end if

    end do

  end do

  do i = 1, nface
    if ( 0 < hvl(i) ) then
      if ( fvl(loc,hvl(i)) <= 0 ) then
        hvl(i) = 0
      end if
    end if
  end do

  f = 0

  do i = 1, nface

    j = hvl(i)

    if ( 0 < j ) then

      f = f + 1

      if ( f < i ) then

        hvl(f) = j

        do

          fvl(facn,j) = f
          j = fvl(succ,j)

          if ( j == hvl(f) ) then
            exit
          end if

        end do

        nrml(1,f) = nrml(1,i)
        nrml(2,f) = nrml(2,i)
        nrml(3,f) = nrml(3,i)

      end if

    end if

  end do

  nface = f

110 continue

  if ( fvl(loc,nvert) <= 0 ) then
    nvert = nvert - 1
    go to 110
  end if

  return
end
