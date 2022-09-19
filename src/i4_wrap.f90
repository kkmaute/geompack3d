function i4_wrap ( ival, ilo, ihi )

!*****************************************************************************80
!
!! I4_WRAP forces an I4 to lie between given limits by wrapping.
!
!  Example:
!
!    ILO = 4, IHI = 8
!
!    I  I4_WRAP
!
!    -2     8
!    -1     4
!     0     5
!     1     6
!     2     7
!     3     8
!     4     4
!     5     5
!     6     6
!     7     7
!     8     8
!     9     4
!    10     5
!    11     6
!    12     7
!    13     8
!    14     4
!
!  Modified:
!
!    15 July 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) IVAL, an integer value.
!
!    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds for the integer value.
!
!    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of IVAL.
!
  implicit none

  integer ( kind = 4 ) i4_modp
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) ihi
  integer ( kind = 4 ) ilo
  integer ( kind = 4 ) ival
  integer ( kind = 4 ) wide

  wide = ihi + 1 - ilo

  if ( wide == 0 ) then
    i4_wrap = ilo
  else
    i4_wrap = ilo + i4_modp ( ival-ilo, wide )
  end if

  return
end
