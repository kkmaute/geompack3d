subroutine updatk ( k, ind, d, e, i, n, p, front, back, fc, ht, ierr )

!*****************************************************************************80
!
!! UPDATK updates a record in FC after a local transformation.
!
!  Discussion:
!
!    This routine updates a record in FC due to a local transformation.
!
!    A simplex with face IND(1:K) and vertex D has D changed to E.
!    Add face IND(1:K) to queue if it is interior face, not yet in
!    queue, and its largest index isn't I.
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
!    Input, integer ( kind = 4 ) K, the number of vertices in a face.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), the vertex indices of face (in any order).
!
!    Input, integer ( kind = 4 ) D, E, the (K+1)st vertex indices of old and new simplices.
!
!    Input, integer ( kind = 4 ) I, the vertex index determining whether face put on queue.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FRONT, BACK, the front and back pointers of queue.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:*), the array of face records;
!    see routine DTRISK.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) p

  integer ( kind = 4 ) back
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) front
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pos

  ierr = 0
  pos = htsrck(k,ind,n,p,fc,ht)

  if ( pos <= 0 ) then
    ierr = 400
    return
  end if

  if ( fc(k+1,pos) == d ) then
    fc(k+1,pos) = e
  else
    fc(k+2,pos) = e
  end if

  if ( fc(k+4,pos) == -1 .and. fc(k,pos) /= i .and. 0 < fc(k+2,pos) ) then

    fc(k+4,pos) = 0

    if ( front == 0 ) then
      front = pos
    else
      fc(k+4,back) = pos
    end if

    back = pos

  end if

  return
end
