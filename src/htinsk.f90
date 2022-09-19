subroutine htinsk ( k, pos, ind, d, e, n, p, fc, ht )

!*****************************************************************************80
!
!! HTINSK inserts a record into the hash table.
!
!  Discussion:
!
!    This routine inserts record FC(1:K+4,POS) containing IND(1:K),D,E,
!    HTLINK,-1 into hash table HT.
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
!    Input, integer ( kind = 4 ) K, number of vertices in a face.
!
!    Input, integer ( kind = 4 ) POS, position of FC array.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), vertex indices of face.  On output,
!    sorted into nondecreasing order.
!
!    Input, integer ( kind = 4 ) D, E, fields K+1, K+2 of FC record (or column).
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:*), array of face records; see routine
!    DTRISK.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) p

  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pos

  call orderk ( k, ind )

  h = ind(1)
  do i = 2, k
    h = mod ( h * n + ind(i), p )
  end do

  fc(1:k,pos) = ind(1:k)
  fc(k+1,pos) = d
  fc(k+2,pos) = e
  fc(k+3,pos) = ht(h)
  fc(k+4,pos) = -1
  ht(h) = pos

  return
end
