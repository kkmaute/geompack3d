subroutine dsftdw ( l, u, k, lda, a, map )

!*****************************************************************************80
!
!! DSFTDW does one step of the heap sort algorithm for double precision data.
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
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) L, U, the lower and upper index of part of heap.
!
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), see routine DHPSRT.
!
!    Input/output, integer ( kind = 4 ) MAP(1:*), see routine DHPSRT.
!
  implicit none

  integer ( kind = 4 ) lda

  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) map(*)
  integer ( kind = 4 ) u

  real    ( kind = 8 ) a(lda,*)
  logical              dless
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) t

  i = l
  j = 2 * i
  t = map(i)

  do

    if ( u < j ) then
      exit
    end if

    if ( j < u ) then
      if ( dless ( k, a(1,map(j)), a(1,map(j+1)) ) ) then
        j = j + 1
      end if
    end if

    if ( dless ( k, a(1,map(j)), a(1,t)) ) then
      exit
    end if

    map(i) = map(j)
    i = j
    j = 2 * i

  end do

  map(i) = t

  return
end
