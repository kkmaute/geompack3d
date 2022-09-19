subroutine rotiar ( n, arr, shift )

!*****************************************************************************80
!
!! ROTIAR rotates elements of an integer array.
!
!  Discussion:
!
!    This routine rotates the elements of an integer array.
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
!   Input, integer ( kind = 4 ) N, the number of elements of array.
!
!   Input/output, integer ( kind = 4 ) ARR(0:N-1), the integer array, which has been
!   shifted on output.
!
!   Input, integer ( kind = 4 ) SHIFT, the amount of (left) shift or rotation; ARR(SHIFT)
!   on input becomes ARR(0) on output.
!
  implicit none

  integer ( kind = 4 ) n

  integer ( kind = 4 ) a
  integer ( kind = 4 ) arr(0:n-1)
  integer ( kind = 4 ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) r
  integer ( kind = 4 ) sh
  integer ( kind = 4 ) shift
  integer ( kind = 4 ) t

  sh = mod ( shift, n )

  if ( sh < 0 ) then
    sh = sh + n
  end if

  if ( sh == 0 ) then
    return
  end if

  a = n
  b = sh

  do

    r = mod ( a, b )
    a = b
    b = r

    if ( r <= 0 ) then
      exit
    end if

  end do

  m = n / a - 1

  do i = 0, a-1
    t = arr(i)
    k = i
    do j = 1, m
      l = k + sh
      if ( n <= l ) then
        l = l - n
      end if
      arr(k) = arr(l)
      k = l
    end do
    arr(k) = t
  end do

  return
end
