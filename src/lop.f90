subroutine lop ( itr, ind, mxedg, top, ldv, vcl, til, tedg, sptr )

!*****************************************************************************80
!
!! LOP applies the local optimization procedure to two triangles.
!
!  Discussion:
!
!    This routine applies a local optimization procedure to two triangles
!    indicated by ITR(1) and ITR(2).  This may result in swapping
!    the diagonal edge of the quadrilateral.
!
!  Modified:
!
!    12 July 1999
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) ITR(1:2), the indices of triangles for LOP.
!
!    Input, integer ( kind = 4 ) IND(1:2), indices indicating common edge of triangles.
!
!    Input, integer ( kind = 4 ) MXEDG, the maximum index of edge to be considered for LOP.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of SPTR indicating top of stack.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TEDG(1:3,1:*), the triangle edge indices;
!    see routine CVDTRI.
!
!    Input/output, integer ( kind = 4 ) SPTR(1:*), stack pointers; see routine CVDTRI.
!
  implicit none

  integer ( kind = 4 ) ldv

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) c
  integer ( kind = 4 ) d
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i4_wrap
  integer ( kind = 4 ) iedg
  integer ( kind = 4 ) in
  integer ( kind = 4 ) ind(2)
  integer ( kind = 4 ) ind1m1
  integer ( kind = 4 ) ind1p1
  integer ( kind = 4 ) ind2m1
  integer ( kind = 4 ) ind2p1
  integer ( kind = 4 ) itr(2)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) mxedg
  integer ( kind = 4 ) top
  integer ( kind = 4 ) sptr(*)
  integer ( kind = 4 ) tedg(3,*)
  integer ( kind = 4 ) til(3,*)
  real    ( kind = 8 ) vcl(ldv,*)
!
!  Common edge is BC, other two vertices are A and D.
!
  iedg = tedg(ind(1),itr(1))
  sptr(iedg) = -1

  ind1m1 = i4_wrap ( ind(1) - 1, 1, 3 )
  ind1p1 = i4_wrap ( ind(1) + 1, 1, 3 )
  ind2m1 = i4_wrap ( ind(2) - 1, 1, 3 )
  ind2p1 = i4_wrap ( ind(2) + 1, 1, 3 )

  b = til(ind(1),itr(1))
  c = til(ind1p1,itr(1))
  a = til(ind1m1,itr(1))
  d = til(ind2m1,itr(2))

  in = diaedg ( vcl(1,d), vcl(2,d), vcl(1,c), vcl(2,c), vcl(1,a), vcl(2,a), &
    vcl(1,b), vcl(2,b) )

  if ( in == 1 ) then
!
!  Check if four edges of quadrilateral should be put on LOP
!  stack, and swap edge BC for AD.
!
   i = tedg(ind1m1,itr(1))

   do j = 1, 4

      if ( j == 2 ) then
        i = tedg(ind1p1,itr(1))
      else if ( j == 3 ) then
        i = tedg(ind2m1,itr(2))
      else if ( j == 4 ) then
        i = tedg(ind2p1,itr(2))
      end if

      if ( i <= mxedg ) then
        if ( sptr(i) == -1 ) then
          sptr(i) = top
          top = i
        end if
      end if

    end do

    til(ind1p1,itr(1)) = d
    til(ind2p1,itr(2)) = a
    tedg(ind(1),itr(1)) = tedg(ind2p1,itr(2))
    tedg(ind(2),itr(2)) = tedg(ind1p1,itr(1))
    tedg(ind1p1,itr(1)) = iedg
    tedg(ind2p1,itr(2)) = iedg

  end if

  return
end
