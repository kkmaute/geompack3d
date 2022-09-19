function htsrck ( k, ind, n, p, fc, ht )

!*****************************************************************************80
!
!! HTSRCK searches for a record in the hash table.
!
!  Discussion:
!
!    This routine searches for record FC(1:K+4,POS) containing key IND(1:K)
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
!    Input, integer ( kind = 4 ) K, number of vertices in a face.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), vertex indices of face.  On output, these
!    have been sorted into nondecreasing order.
!
!    Input, integer ( kind = 4 ) N, upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, size of hash table.
!
!    Input, integer ( kind = 4 ) FC(1:K+4,1:*), array of face records; see routine DTRISK.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), hash table using direct chaining.
!
!    Output, integer ( kind = 4 ) HTSRCK, position of FC record with key IND(1:K) if found,
!    or 0 if not found.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) kp3
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pos

  kp3 = k + 3
  call orderk ( k, ind )

  h = ind(1)
  do i = 2, k
    h = mod ( h * n + ind(i), p )
  end do

  pos = ht(h)

20 continue

  if ( pos /= 0 ) then

    i = 1

30  continue

    if ( fc(i,pos) /= ind(i) ) then
      pos = fc(kp3,pos)
      go to 20
    end if

    i = i + 1

    if ( i <= k ) then
      go to 30
    end if

  end if

  htsrck = pos

  return
end
