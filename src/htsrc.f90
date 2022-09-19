function htsrc ( a, b, c, n, p, fc, ht )

!*****************************************************************************80
!
!! HTSRC searches for a record in the hash table.
!
!  Discussion:
!
!    This routine searches for record FC(1:7,IND) containing key A,B,C
!    in hash table HT.
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
!    Input, integer ( kind = 4 ) A,B,C, first 3 fields of FC record (in any order).
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input, integer ( kind = 4 ) FC(1:7,1:*), array of face records; see routine DTRIS3.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) HTSRC, index of FC record with key A,B,C if found,
!    or 0 if not found.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  integer ( kind = 4 ) aa
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bb
  integer ( kind = 4 ) c
  integer ( kind = 4 ) cc
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrc
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n

  aa = a
  bb = b
  cc = c
  call order3 ( aa, bb, cc )
  k = mod ( aa * n + bb, p )
  k = mod ( k * n + cc, p )
  ind = ht(k)

  do

    if ( ind == 0 ) then
      exit
    end if

    if ( fc(1,ind) == aa .and. fc(2,ind) == bb .and. fc(3,ind) == cc ) then
      exit
    end if

    ind = fc(6,ind)

  end do

  htsrc = ind

  return
end
