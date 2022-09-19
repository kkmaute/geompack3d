subroutine rmcled ( nface, nvert, hvl, fvl )

!*****************************************************************************80
!
!! RMCLED removes collinear adjacent convex polyhedron edges from the database.
!
!  Discussion:
!
!    This routine removes collinear adjacent edges from a convex polyhedron
!    data structure (assuming no coplanar adjacent faces).
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
!    Input, integer ( kind = 4 ) NFACE, the number of faces in convex polyhedron.
!
!    Input/output, integer ( kind = 4 ) NVERT, the size of FVL; must be greater than
!    or equal to twice the number of edges of polyhedron.
!
!    Input/output, integer ( kind = 4 ) HVL(1:NFACE), head vertex list.
!
!    Input/output, integer ( kind = 4 ) FVL(1:5,1:NVERT), face vertex list; see routine
!    DSCPH; may contain unused columns, indicated by LOC values <= 0.
!
  implicit none

  integer ( kind = 4 ) nface
  integer ( kind = 4 ) nvert

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ), parameter :: edgv = 5
  integer ( kind = 4 ) f
  integer ( kind = 4 ), parameter :: facn = 2
  integer ( kind = 4 ) fvl(5,nvert)
  integer ( kind = 4 ) g
  integer ( kind = 4 ) hvl(nface)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ), parameter :: succ = 3

  do f = 1, nface

    j = hvl(f)

    do

      i = fvl(pred,j)
      a = fvl(edgv,j)
      b = fvl(edgv,i)

      if ( fvl(facn,a) /= fvl(facn,b) ) then
        exit
      end if

      j = fvl(succ,j)

    end do

    hvl(f) = j
    i = j
    b = fvl(edgv,i)

30  continue

    j = fvl(succ,i)
    a = fvl(edgv,j)

    if ( fvl(facn,a) == fvl(facn,b) ) then

      g = fvl(facn,b)
      c = fvl(succ,b)

      do

        fvl(loc,j) = 0
        fvl(loc,b) = 0
        j = fvl(succ,j)
        b = a
        a = fvl(edgv,j)

        if ( fvl(facn,a) /= g ) then
          exit
        end if

      end do

      fvl(succ,i) = j
      fvl(pred,j) = i
      fvl(succ,b) = c
      fvl(pred,c) = b
      fvl(edgv,i) = b
      fvl(edgv,b) = i

    end if

    i = j
    b = a

    if ( i /= hvl(f) ) then
      go to 30
    end if

  end do

  do while ( fvl(loc,nvert) <= 0 )
    nvert = nvert - 1
  end do

  return
end
