subroutine mtredg ( utype, i1, i2, i3, ibndry, nt, til, tedg )

!*****************************************************************************80
!
!! MTREDG sets fields for a triangle as needed by routine TMERGE.
!
!  Discussion:
!
!    This routine sets fields for triangle as needed by routine TMERGE.
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
!    Input, logical UTYPE, is TRUE iff triangle contains two 'U' vertices.
!
!    Input, integer ( kind = 4 ) I1, I2, I3, the indices of 3 triangle vertices in VCL;
!    the first two indices also belong to the next merge edge.
!
!    Input, integer ( kind = 4 ) IBNDRY, the index of boundary edge for TEDG.
!
!    Input/output, integer ( kind = 4 ) NT, the number of entries in TIL, TEDG so far.
!
!    Input/output, integer ( kind = 4 ) TIL(1:NT), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TEDG(1:NT), the triangle edge indices; see
!    routine TMERGE.
!
  implicit none

  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) ibndry
  integer ( kind = 4 ) nt
  integer ( kind = 4 ) tedg(3,*)
  integer ( kind = 4 ) til(3,*)
  logical              utype

  nt = nt + 1
  til(1,nt) = i1
  til(2,nt) = i2
  til(3,nt) = i3
  tedg(1,nt) = nt

  if ( utype ) then
    tedg(2,nt) = nt - 1
    tedg(3,nt) = ibndry
  else
    tedg(2,nt) = ibndry
    tedg(3,nt) = nt - 1
  end if

  return
end
