subroutine ihpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! IHPSRT sorts a list of integer points in KD.
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    integer ( kind = 4 ) points so that the points are in lexicographic
!    increasing order.
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
!    Input, integer ( kind = 4 ) N, number of points.
!
!    Input, integer ( kind = 4 ) LDA, leading dimension of array A in calling routine;
!    should be greater than or equal to K.
!
!    Input, integer ( kind = 4 ) A(1:K,1:*), array of at least N K-D integer points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, the points of A with indices
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, elements
!    are permuted so that A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N))
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  integer ( kind = 4 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call isftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call isftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
