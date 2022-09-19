subroutine randpt ( k, n, seed, axis, nptav, scale, trans, lda, a )

!*****************************************************************************80
!
!! RANDPT generates N random points in KD.
!
!  Discussion:
!
!    This routine generates N random K-dimensional points from the uniform
!    distribution.
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
!    Input, integer ( kind = 4 ) K, the dimension of points.
!
!    Input, integer ( kind = 4 ) N, the number of random points.
!
!    Input/output, integer ( kind = 4 ) SEED, the seed for pseudo random number generator.
!
!    Input, integer ( kind = 4 ) AXIS, NPTAV; if AXIS < 1 or K < AXIS, then uniform
!    random points are generated; if 1 <= AXIS <= K then an average of NPTAV
!    uniform random points are generated with the same AXIS
!    coordinate on about N/NPTAV random parallel hyperplanes.
!
!    Input, real ( kind = 8 ) SCALE(1:K), TRANS(1:K), scale and translation
!    factors for coordinates 1 to K; Ith coordinate of random point is
!    R*SCALE(I) + TRANS(I) where 0 < R < 1.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of array A in calling
!    routine; should be greater than or equal to K.
!
!    Output, real ( kind = 8 ) A(1:K,1:N), an array of N uniform random
!    K-D points.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) n

  real    ( kind = 8 ) a(lda,n)
  integer ( kind = 4 ) axis
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nptav
  real    ( kind = 8 ) r
  real    ( kind = 8 ) scale(k)
  integer ( kind = 4 ) seed
  real    ( kind = 8 ) trans(k)
  real    ( kind = 8 ) urand

  if ( axis < 1 .or. k < axis ) then

    do j = 1, n
      do i = 1, k
        a(i,j) = urand ( seed ) * scale(i) + trans(i)
      end do
    end do

  else

    m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
    r = urand ( seed ) * scale(axis) + trans(axis)

    do j = 1, n
      do i = 1, k
        if ( i == axis ) then
          a(i,j) = r
        else
          a(i,j) = urand ( seed ) * scale(i) + trans(i)
        end if
      end do

      m = m - 1

      if ( m <= 0 ) then
        m = int ( urand ( seed ) * 2.0D+00 * nptav + 0.5D+00 )
        r = urand ( seed ) * scale(axis) + trans(axis)
      end if

    end do

  end if

  return
end
