subroutine edght ( a, b, v, n, htsiz, maxedg, hdfree, last, ht, edge, w, ierr )

!*****************************************************************************80
!
!! EDGHT searches a hash table for an edge record.
!
!  Discussion:
!
!    This routine searches in the hash table HT for the record in EDGE
!    containing the key (A,B).
!
!    Before first call to this routine, HDFREE, LAST, and
!    entries of HT should be set to 0.
!
!  Modified:
!
!    17 August 2005
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
!    Input, integer ( kind = 4 ) A, B, vertex indices, greater than 0, of edge (also
!    key of hash table).
!
!    Input, integer ( kind = 4 ) V, value associated with edge.
!
!    Input, integer ( kind = 4 ) N, upper bound on A, B.
!
!    Input, integer ( kind = 4 ) HTSIZ, size of hash table HT.
!
!    Input, integer ( kind = 4 ) MAXEDG, maximum size available for EDGE array.
!
!    Input/output, integer ( kind = 4 ) HDFREE, head pointer to linked list of free
!    entries of EDGE array due to deletions.
!
!    Input/output, integer ( kind = 4 ) LAST, index of last entry used in EDGE array.
!
!    Input/output, integer ( kind = 4 ) HT(0:HTSIZ-1), hash table of head pointers
!    (direct chaining with ordered lists is used).  If key with A,B is
!    found then this record is deleted from hash table, else record is
!    inserted in hash table
!
!    Input/output, integer ( kind = 4 ) EDGE(1:4,1:MAXEDG), entries of hash table records;
!    EDGE(1,I) = MIN ( A, B ); EDGE(2,I) = MAX ( A, B );
!    EDGE(3,I) = V; EDGE(4,I) = link.
!    If key with A,B is found then this record is deleted
!    from hash table, else record is inserted in hash table
!
!    Output, integer ( kind = 4 ) W, EDGE(3,INDEX), where INDEX is index of record,
!    if found; else 0.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) htsiz
  integer ( kind = 4 ) maxedg

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) bptr
  integer ( kind = 4 ) edge(4,maxedg)
  integer ( kind = 4 ) hdfree
  integer ( kind = 4 ) ht(0:htsiz-1)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) k
  integer ( kind = 4 ) last
  integer ( kind = 4 ) n
  integer ( kind = 4 ) newp
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) v
  integer ( kind = 4 ) w

  ierr = 0

  if ( a < b ) then
    aa = a
    bb = b
  else
    aa = b
    bb = a
  end if

  k = mod ( aa * n + bb, htsiz )
  bptr = -1
  ptr = ht(k)

  do while ( ptr /= 0 )

    if ( aa < edge(1,ptr) ) then

      exit

    else if ( edge(1,ptr) == aa ) then

      if ( bb < edge(2,ptr) ) then

        exit

      else if ( edge(2,ptr) == bb ) then

        if ( bptr == -1 ) then
          ht(k) = edge(4,ptr)
        else
          edge(4,bptr) = edge(4,ptr)
        end if

        edge(4,ptr) = hdfree
        hdfree = ptr
        w = edge(3,ptr)
        return

      end if

    end if

    bptr = ptr
    ptr = edge(4,ptr)

  end do

  if ( 0 < hdfree ) then
    newp = hdfree
    hdfree = edge(4,hdfree)
  else
    last = last + 1
    newp = last
    if ( maxedg < last ) then
      ierr = 1
      return
    end if
  end if

  if ( bptr == -1 ) then
    ht(k) = newp
  else
    edge(4,bptr) = newp
  end if

  edge(1,newp) = aa
  edge(2,newp) = bb
  edge(3,newp) = v
  edge(4,newp) = ptr
  w = 0

  return
end
