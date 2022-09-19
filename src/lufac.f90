subroutine lufac ( a, lda, n, tol, ipvt, singlr )

!*****************************************************************************80
!
!! LUFAC factors a matrix.
!
!  Discussion:
!
!    This routine obtains the LU factorization of a matrix A by applying
!    Gaussian elimination with partial pivoting.
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
!    Input/output, real ( kind = 8 ) A(LDA,N).  On input, the N by N matrix
!    to be factored.  On output, the upper triangular matrix U and multipliers
!    of unit lower triangular matrix L (if matrix A is nonsingular).
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling routine.
!
!    Input, integer ( kind = 4 ) N, the order of matrix A.
!
!    Input, real ( kind = 8 ) TOL, the relative tolerance for detecting
!    singularity of A.
!
!    Output, integer ( kind = 4 ) IPVT(1:N-1), the pivot indices.
!
!    Output, logical SINGLR, TRUE iff matrix is singular; this occurs when the
!    magnitude of a pivot element is <= TOL * MAX ( |A(I,J)| ).
!
  implicit none

  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ipvt(n-1)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) m
  logical              singlr
  real    ( kind = 8 ) t
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs

  if ( n < 1 ) then
    return
  end if

  singlr = .true.

  t = maxval ( abs ( a(1:n,1:n) ) )

  tolabs = tol * t

  do k = 1, n-1

    kp1 = k + 1
    m = k

    do i = k+1, n
      if ( abs ( a(m,k) ) < abs ( a(i,k) ) ) then
        m = i
      end if
    end do

    ipvt(k) = m

    t = a(m,k)
    a(m,k) = a(k,k)
    a(k,k) = t

    if ( abs ( t ) <= tolabs ) then
      return
    end if

    a(kp1:n,k) = a(kp1:n,k) / t

    do j = kp1, n
      t = a(m,j)
      a(m,j) = a(k,j)
      a(k,j) = t
      a(kp1:n,j) = a(kp1:n,j) - a(kp1:n,k) * t
    end do

  end do

  if ( tolabs < abs ( a(n,n) ) ) then
    singlr = .false.
  end if

  return
end
