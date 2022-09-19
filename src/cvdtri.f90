subroutine cvdtri ( inter, ldv, nt, vcl, til, tedg, sptr, ierror )

!*****************************************************************************80
!
!! CVDTRI converts boundary triangles to Delaunay triangles.
!
!  Discussion:
!
!    This routine converts triangles in a strip near the boundary of
!    the polygon or inside the polygon to Delaunay triangles.
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
!  Reference:
!
!    Barry Joe,
!    GEOMPACK - a software package for the generation of meshes
!      using geometric algorithms,
!    Advances in Engineering Software,
!    Volume 13, pages 325-331, 1991.
!
!  Parameters:
!
!    Input, logical INTER, is TRUE if and only if there is at least
!    one interior mesh vertex.
!
!    Input, integer ( kind = 4 ) LDV, the leading dimension of VCL in calling routine.
!
!    Input, integer ( kind = 4 ) NT, the number of triangles in strip or polygon.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:NT), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TEDG(1:3,1:NT), TEDG(J,I) refers to edge with
!    vertices TIL(J:J+1,I) and contains index of merge edge or
!    greater than NT for edge of chains.
!
!    Workspace, integer SPTR(1:NT); SPTR(I) = -1 if merge edge I is not in
!    LOP stack, else greater than or equal to 0 and pointer (index of SPTR)
!    to next edge in stack (0 indicates bottom of stack).
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) ldv
  integer ( kind = 4 ) nt

  integer ( kind = 4 ) e
  integer ( kind = 4 ) ierror
  integer ( kind = 4 ) ind(2)
  logical              inter
  integer  ( kind = 4 ) itr(2)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) mxtr
  logical              sflag
  integer ( kind = 4 ) sptr(nt)
  integer ( kind = 4 ) tedg(3,nt)
  integer ( kind = 4 ) til(3,nt)
  integer ( kind = 4 ) top
  real    ( kind = 8 ) vcl(ldv,*)

  ierror = 0
  sflag = .true.
  sptr(1:nt) = -1

  do k = 1, nt

    mxtr = k + 1

    if ( k == nt ) then
      if ( .not. inter ) then
        return
      end if
      mxtr = nt
      sflag = .false.
    end if

    top = k
    sptr(k) = 0

    do

      e = top
      top = sptr(e)

      call fndtri ( e, mxtr, sflag, tedg, itr, ind, ierror )

      if ( ierror /= 0 ) then
        return
      end if

      call lop ( itr, ind, k, top, ldv, vcl, til, tedg, sptr )

      if ( top <= 0 ) then
        exit
      end if

    end do

  end do

  return
end
