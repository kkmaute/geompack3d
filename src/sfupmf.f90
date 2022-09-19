subroutine sfupmf ( r, psi, indp, loch )

!*****************************************************************************80
!
!! SFUPMF sifts up a heap.
!
!  Discussion:
!
!    This routine sifts PSI(INDP(R)) up a heap which has maximum PSI value
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
!    Input, integer ( kind = 4 ) R, the element of heap to be sifted up.
!
!    Input, real ( kind = 8 ) PSI(1:*), the key values for heap.
!
!    Input/output, integer ( kind = 4 ) INDP(1:R), the indices of PSI which are maintained
!    in heap.
!
!    Input/output, integer ( kind = 4 ) LOCH(1:*), the location of indices in heap
!    (inverse of INDP).
!
  implicit none

  integer ( kind = 4 ) r

  integer ( kind = 4 ) i
  integer ( kind = 4 ) indp(r)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) loch(*)
  real    ( kind = 8 ) psi(*)
  real    ( kind = 8 ) t

  i = r
  j = int ( i / 2 )
  k = indp(i)
  t = psi(k)

  do

    if ( i <= 1 ) then
      exit
    end if

    if ( t <= psi(indp(j)) ) then
      exit
    end if

    indp(i) = indp(j)
    loch(indp(i)) = i
    i = j
    j = int ( i / 2 )

  end do

  indp(i) = k
  loch(k) = i

  return
end
