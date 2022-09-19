subroutine dhpsrt ( k, n, lda, a, map )

!*****************************************************************************80
!
!! DHPSRT sorts a list of double precision points in KD.
!
!  Discussion:
!
!    This routine uses heapsort to obtain the permutation of N K-dimensional
!    double precision points so that the points are in lexicographic
!    increasing order.
!
!  Modified:
!
!    31 August 2005
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
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of points.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling
!    routine; K <= LDA.
!
!    Input, real ( kind = 8 ) A(1:K,1:*), array of points.
!
!    Input/output, integer ( kind = 4 ) MAP(1:N).  On input, he points of A with indices
!    MAP(1), MAP(2), ..., MAP(N) are to be sorted.  On output, the elements
!    are permuted so that A(*,MAP(1)) <= A(*,MAP(2)) <= ... <= A(*,MAP(N)).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) map(n)
  integer ( kind = 4 ) t

  do i = n/2, 1, -1
    call dsftdw ( i, n, k, lda, a, map )
  end do

  do i = n, 2, -1
    t = map(1)
    map(1) = map(i)
    map(i) = t
    call dsftdw ( 1, i-1, k, lda, a, map )
  end do

  return
end
