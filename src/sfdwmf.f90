subroutine sfdwmf ( l, r, psi, indp, loch )

!*****************************************************************************80
!
!! SFDWMF sifts down a heap.
!
!  Discussion:
!
!    This routine sifts PSI(INDP(L)) down a heap which has maximum PSI value
!    at the root of the heap and is maintained by pointers in INDP.
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
!    Input, integer ( kind = 4 ) L, element of heap to be sifted down.
!
!    Input, integer ( kind = 4 ) R, upper bound of heap.
!
!    Input, real ( kind = 8 ) PSI(1:*), key values for heap.
!
!    Input/output, integer ( kind = 4 ) INDP(1:R), indices of PSI which are maintained
!    in heap.
!
!    Input/output, integer ( kind = 4 ) LOCH(1:*), location of indices in heap (inverse
!    of INDP).
!
  implicit none

  integer ( kind = 4 ) r

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(r)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) loch(*)
  real    ( kind = 8 ) psi(*)
  real    ( kind = 8 ) t

  i = l
  j = 2 * i
  k = indp(i)
  t = psi(k)

  do

    if ( r < j ) then
      exit
    end if

    if ( j < r ) then
      if ( psi(indp(j)) < psi(indp(j+1)) ) then
        j = j + 1
      end if
    end if

    if ( psi(indp(j)) <= t ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = 2 * i

  end do

  indp(i) = k
  loch(k) = i

  return
end
