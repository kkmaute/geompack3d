subroutine htdelk ( k, pos, n, p, fc, ht )

!*****************************************************************************80
!
!! HTDELK deletes a record from the hash table.
!
!  Discussion:
!
!    This routine deletes record FC(1:K+4,POS) from hash table HT.
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
!    Input, integer ( kind = 4 ) POS, the position of FC array.
!
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:*), the array of face records;
!    see routine DTRISK.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) h
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) kp3
  integer ( kind = 4 ) n
  integer ( kind = 4 ) pos
  integer ( kind = 4 ) ptr

  kp3 = k + 3
  h = fc(1,pos)

  do i = 2, k
    h = mod ( h*n + fc(i,pos), p )
  end do

  ptr = ht(h)

  if ( ptr == pos ) then

    ht(h) = fc(kp3,pos)

  else

    do while ( fc(kp3,ptr) /= pos )
      ptr = fc(kp3,ptr)
    end do

    fc(kp3,ptr) = fc(kp3,pos)

  end if

  return
end
