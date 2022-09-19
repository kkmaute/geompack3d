subroutine orderk ( k, ind )

!*****************************************************************************80
!
!! ORDERK reorders K elements of an array in nondecreasing order.
!
!  Discussion:
!
!    This routine reorders K elements of array IND in nondecreasing order.
!
!    It is assumed that K is small, say <= 15, so that insertion sort
!    is used. If K is larger, a faster sort such as heapsort should
!    be used.
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
!    Input, integer ( kind = 4 ) K, the size of array IND.
!
!    Input/output, integer ( kind = 4 ) IND(1:K), an array, which is sorted on output.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) s
  integer ( kind = 4 ) t

  do i = 2, k

    t = ind(i)
    j = i

    do

      s = ind(j-1)

      if ( s <= t ) then
        exit
      end if

      ind(j) = s
      j = j - 1

      if ( j <= 1 ) then
        exit
      end if

    end do

    ind(j) = t

  end do

  return
end
