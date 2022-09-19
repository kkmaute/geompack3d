subroutine htdel ( ind, n, p, fc, ht )

!*****************************************************************************80
!
!! HTDEL deletes a record from the hash table.
!
!  Discussion:
!
!    This routine deletes record FC(1:7,IND) from the hash table HT.
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
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records;
!    see routine DTRIS3.  On output, one link in FC is updated.
!
!    Input/output, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!    On output, one link in HT is updated.
!
  implicit none

  integer ( kind = 4 ) p

  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) ind
  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ptr

  k = mod ( fc(1,ind) * n + fc(2,ind), p )
  k = mod ( k * n + fc(3,ind), p )
  ptr = ht(k)

  if ( ptr == ind ) then

    ht(k) = fc(6,ind)

  else

    do while ( fc(6,ptr) /= ind )
      ptr = fc(6,ptr)
    end do
    fc(6,ptr) = fc(6,ind)

  end if

  return
end
