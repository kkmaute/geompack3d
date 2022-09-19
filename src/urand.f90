function urand ( iy )

!*****************************************************************************80
!
!! URAND is a uniform random number generator.
!
!  Discussion:
!
!    This routine is a uniform random number generator based on theory and
!    suggestions given in D. E. Knuth (1969), Vol. 2. The integer IY
!    should be initialized to an arbitrary integer prior to the first
!    call to URAND. The calling program should not alter the value of
!    IY between subsequent calls to URAND. Values of URAND will be
!    returned in the interval (0,1).
!
!  Reference:
!
!    George Forsythe, Michael Malcolm, Cleve Moler,
!    Computer Methods for Mathematical Computations,
!    Prentice Hall, 1971, page 246.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) IY, a seed for the generator.
!
!    Output, real ( kind = 8 ) URAND, a pseudorandom value.
!
  implicit none

  real    ( kind = 8 ) halfm
  integer ( kind = 4 ), save :: ia = 0
  integer ( kind = 4 ), save :: ic = 0
  integer ( kind = 4 ) iy
  integer ( kind = 4 ) m
  integer ( kind = 4 ), save :: m2 = 0
  integer ( kind = 4 ), save :: mic = 0
  real    ( kind = 8 ), save :: s = 0.0E+00
  real    ( kind = 8 ) urand

  if ( m2 == 0 ) then
!
!  If first entry, compute machine integer word length
!
    m = 1

    do

      m2 = m
      m = 2 * m2

      if ( m <= m2 ) then
        exit
      end if

    end do

    halfm = m2
!
!  Compute multiplier and increment for linear congruential method
!
    ia = 8 * int ( halfm * atan (1.0D+00 ) / 8.0D+00 ) + 5
    ic = 2 * int ( halfm * ( 0.5D+00 - sqrt ( 3.0D+00 ) / 6.0D+00 ) ) + 1
    mic = ( m2 - ic ) + m2
!
!  S is the scale factor for converting to floating point
!
    s = 0.5D+00 / halfm

  end if
!
!  Compute next random number
!
  iy = iy * ia
!
!  The following statement is for computers which do not allow
!  integer ( kind = 4 ) overflow on addition
!
  if ( mic < iy ) then
    iy = ( iy - m2 ) - m2
  end if

  iy = iy + ic
!
!  The following statement is for computers where the word
!  length for addition is greater than for multiplication
!
  if ( m2 < iy / 2 ) then
    iy = (iy - m2) - m2
  end if
!
!  The following statement is for computers where integer
!  overflow affects the sign bit
!
  if ( iy < 0 ) then
    iy = ( iy + m2 ) + m2
  end if

  urand = real ( iy, kind = 8 ) * s

  return
end
