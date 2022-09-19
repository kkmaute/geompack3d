subroutine lusol ( a, lda, n, ipvt, b )

!*****************************************************************************80
!
!! LUSOL solves a linear system involving a matrix factored by LUFAC.
!
!  Discussion:
!
!    This routine solves a linear system A*X = B given LU factorization of A.
!    It is assumed that A is nonsingular.
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
!    Input, real ( kind = 8 ) A(1:N,1:N), contains factors L, U output
!    from routine LUFAC.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices from routine LUFAC.
!
!    Input/output, real ( kind = 8 ) B(1:N).  On input, the right hand
!    side vector.  On output, the solution vector X
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  real    ( kind = 8 ) b(n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) m
  real    ( kind = 8 ) t
!
!  Forward elimination
!
  do k = 1, n-1
    m = ipvt(k)
    t = b(m)
    b(m) = b(k)
    b(k) = t
    do i = k+1, n
      b(i) = b(i) - a(i,k)*t
    end do

  end do
!
!  Back substitution
!
  do k = n,2,-1
    t = b(k)/a(k,k)
    b(k) = t
    b(1:k-1) = b(1:k-1) - a(1:k-1,k)*t
  end do

  b(1) = b(1) / a(1,1)

  return
end
