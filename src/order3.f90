subroutine order3 ( i, j, k )

!*****************************************************************************80
!
!! ORDER3 reorders 3 integers into ascending order.
!
!  Discussion:
!
!    This routine reorders I, J, K so that I <= J <= K.
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
!    Input/output, integer ( kind = 4 ) I, J, K, on output are sorted into
!    nondecreasing order.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) t

  if ( j < i ) then
    if ( k < j ) then
      call i4_swap ( i, k )
    else if ( k < i ) then
      t = i
      i = j
      j = k
      k = t
    else
      call i4_swap ( i, j )
    end if
  else
    if ( k < i ) then
      t = i
      i = k
      k = j
      j = t
    else if ( k < j ) then
      call i4_swap ( j, k )
    end if
  end if

  return
end
