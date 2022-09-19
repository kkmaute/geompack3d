subroutine stats ( n, x, a, h, nf, xmin, xmax, mean, stdv, freq )

!*****************************************************************************80
!
!! STATS computes statistical measurements for data.
!
!  Discussion:
!
!    This routine computes the following statistical measurements
!    for data X(1:N):
!
!    * minimum,
!    * maximum,
!    * mean,
!    * standard deviation,
!    * relative frequency in intervals
!      X(I) < A+H,
!             A+H <= X(I) < A+2*H, ...,
!                           A+(NF-1)*H <= X(I) < A+NF*H, and
!                                                A+NF*H <= X(I)
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
!    Input, integer ( kind = 4 ) N, the number of data points (measurements).
!
!    Input, real ( kind = 8 ) X(1:N), array of N data points.
!
!    Input, real ( kind = 8 ) A, H, used for determining frequency
!    intervals as above.
!
!    Input, integer ( kind = 4 ) NF, the number of frequency intervals - 1.
!
!    Output, real ( kind = 8 ) XMIN, XMAX, MEAN, STDV, the minimum,
!    maximum, mean, standard deviation.
!
!    Output, real ( kind = 8 ) FREQ(0:NF), the relative frequencies
!    (fraction of 1) as above.
!
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf

  real    ( kind = 8 ) a
  real    ( kind = 8 ) freq(0:nf)
  real    ( kind = 8 ) h
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) stdv
  real    ( kind = 8 ) xmax
  real    ( kind = 8 ) xmin
  real    ( kind = 8 ) x(n)

  freq(0:nf) = 0.0D+00

  xmin = x(1)
  xmax = x(1)

  do i = 1, n
    xmin = min ( xmin, x(i) )
    xmax = max ( xmax, x(i) )
    k = int((x(i) - a)/h)
    k = max ( min ( k, nf ), 0 )
    freq(k) = freq(k) + 1.0D+00
  end do

  mean = sum ( x(1:n) ) / real ( n, kind = 8 )

  stdv = sqrt ( sum ( ( x(1:n) - mean )**2 ) / real ( n-1, kind = 8 ) )

  freq(0:nf) = freq(0:nf) / real ( n, kind = 8 )

  return
end
