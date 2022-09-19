subroutine updatf ( a, b, c, d, e, i, n, p, front, back, fc, ht, ierr )

!*****************************************************************************80
!
!! UPDATF updates a record in FC after a local transformation.
!
!  Discussion:
!
!    This routine updates a record in FC due to a local transformation.
!
!    Tetrahedron ABCD becomes ABCE. Add face ABC to queue if it is
!    interior face, not yet in queue, and its largest index isn't I.
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
!    Input, integer ( kind = 4 ) A, B, C, the first 3 fields of FC record (in any order).
!
!    Input, integer ( kind = 4 ) D, E, the fourth vertex indices of old and new tetrahedrons.
!
!    Input, integer ( kind = 4 ) I, the vertex index determining whether face put on queue.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the front and back pointers of queue.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) back
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) n

  ierr = 0
  ind = htsrc ( a, b, c, n, p, fc, ht )

  if ( ind <= 0 ) then
    ierr = 300
    return
  end if

  if ( fc(4,ind) == d ) then
    fc(4,ind) = e
  else
    fc(5,ind) = e
  end if

  if ( fc(7,ind) == -1 .and. fc(3,ind) /= i .and. 0 < fc(5,ind) ) then

    fc(7,ind) = 0

    if ( front == 0 ) then
      front = ind
    else
      fc(7,back) = ind
    end if

    back = ind
  end if

  return
end
