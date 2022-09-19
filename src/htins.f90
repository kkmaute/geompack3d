subroutine htins ( ind, a, b, c, d, e, n, p, fc, ht )

!*****************************************************************************80
!
!! HTINS inserts a record into the hash table.
!
!  Discussion:
!
!    This routine inserts record FC(1:7,IND) containing A,B,C,D,E,HTLINK,-1
!    into hash table HT.
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
!    Input, integer ( kind = 4 ) IND, the index of FC array.
!
!    Input, integer ( kind = 4 ) A, B, C, D, E, the first 5 fields of FC record (or column).
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records;
!    see routine DTRIS3.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) d
  integer ( kind = 4 ) e
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  fc(1,ind) = aa
  fc(2,ind) = bb
  fc(3,ind) = cc
  fc(4,ind) = d
  fc(5,ind) = e
  fc(6,ind) = ht(k)
  fc(7,ind) = -1
  ht(k) = ind

  return
end
