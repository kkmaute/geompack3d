subroutine isftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! ISFTDW does one step of the heap sort algorithm for integer data.
!
!  Discussion:
!
!    This routine sifts A(*,MAP(L)) down a heap of size U.
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
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) A(1:K,1:*), see routine IHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*), see routine IHPSRT.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) a(lda,*)
  integer ( kind = 4 ) i
  logical              iless
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) u

  i = l
  j = 2 * i
  t = map(i)

  do while ( j <= u )

    if ( j < u ) then
      if ( iless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( iless ( k, a(1,map(j)), a(1,t) ) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2*i

  end do

  map(i) = t

  return
end
