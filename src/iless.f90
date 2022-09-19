function iless ( k, p, q )

!*****************************************************************************80
!
!! ILESS determines the lexicographically lesser of two integer values.
!
!  Discussion:
!
!    This routine determines whether P is lexicographically less than Q.
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
!    Input, integer ( kind = 4 ) K, dimension of points.
!
!    Input, integer ( kind = 4 ) P(1:K), Q(1:K), two K-dimensional integer points.
!
!    Output, logical ILESS, TRUE if P < Q, FALSE otherwise.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) i
  logical              iless
  integer ( kind = 4 ) p(k)
  integer ( kind = 4 ) q(k)

  do i = 1, k

    if ( p(i) == q(i) ) then
      cycle
    end if

    if ( p(i) < q(i) ) then
      iless = .true.
    else
      iless = .false.
    end if

    return

  end do

  iless = .false.

  return
end
