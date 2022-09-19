subroutine fndtri ( iedg, mxtr, sflag, tedg, itr, ind, ierror )

!*****************************************************************************80
!
!! FNDTRI finds two triangles containing a given edge.
!
!  Discussion:
!
!    This routine finds two triangles containing the edge with index IEDG
!    in array TEDG.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IEDG, the index of edge to be searched in TEDG.
!
!    Input, integer ( kind = 4 ) MXTR, the maximum index of triangle to be searched in TEDG.
!
!    Input, logical SFLAG, is TRUE if and only if the second triangle is to be
!    searched from end of array.
!
!    Input, integer ( kind = 4 ) TEDG(1:3,1:MXTR), triangle edge indices; see routine CVDTRI.
!
!    Output, integer ( kind = 4 ) ITR(1:2), IND(1:2), indices such that IEDG =
!    TEDG(IND(1),ITR(1)) = TEDG(IND(2),ITR(2)).
!
!    Output, integer ( kind = 4 ) IERROR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) mxtr

  integer ( kind = 4 ) i
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(2)
  integer ( kind = 4 ) itr(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical              sflag
  integer ( kind = 4 ) tedg(3,mxtr)
!
!  Search from end of array TEDG.
!
  ierror = 0
  k = 1
  j = 1
  i = mxtr

10 continue

  do

    if ( tedg(j,i) == iedg ) then
      exit
    end if

    j = j + 1

    if ( 3 < j ) then
      j = 1
      i = i - 1
      if ( i <= 0 ) then
        ierror = 231
        return
      end if
    end if

  end do

  itr(k) = i
  ind(k) = j

  if ( k == 2 ) then
    return
  end if

  k = 2

  if ( sflag ) then

    j = 1
    i = i - 1

    if ( i <= 0 ) then
      ierror = 231
      return
    end if

    go to 10

  end if
!
!  Search from beginning of array TEDG for second triangle.
!
  j = 1
  i = 1
   20 continue

  if ( itr(1) <= i ) then
    ierror = 231
    return
  end if

   30 continue

  if ( tedg(j,i) /= iedg ) then
    j = j + 1
    if ( 3 < j ) then
      j = 1
      i = i + 1
      go to 20
    else
      go to 30
    end if
  end if

  itr(2) = i
  ind(2) = j

  return
end
