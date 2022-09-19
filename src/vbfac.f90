subroutine vbfac ( pt, ctr, vcl, vm, bf, fc, topv, topnv )

!*****************************************************************************80
!
!! VBFAC determines the boundary faces of a 3D triangulation.
!
!  Discussion:
!
!    This routine determines boundary faces of a 3D triangulation visible
!    from point PT, given a starting visible boundary face.
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
!    Input, real ( kind = 8 ) PT(1:3), the 3D point.
!
!    Input, real ( kind = 8 ) CTR(1:3), the 3D point in interior of
!    triangulation.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL being triangulated.
!
!    Input, integer ( kind = 4 ) BF(1:3,1:*), the array of boundary face records; see DTRIS3.
!
!    Input/output, integer ( kind = 4 ) FC(1:7,1:*), the array of face records; see routine
!    DTRIS3; row 7 is used for links of 3 stacks in this routine.  On output,
!    FC(7,*) has been updated, so that only stack of visible boundary
!    faces remains.
!
!    Input/output, integer ( kind = 4 ) TOPV.  On input, index of FC of visible boundary
!    face.  On output, index of top of stack of visible boundary faces.
!
!    Input, integer ( kind = 4 ) TOPNV, the index of top of stack of boundary faces
!    already found to be not visible from PT, or 0 for empty stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) bf(3,*)
  real    ( kind = 8 ) ctr(3)
  integer ( kind = 4 ) fc(7,*)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) op
  integer ( kind = 4 ) opside
  real    ( kind = 8 ) pt(3)
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topn
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) topt
  integer ( kind = 4 ) topv
  integer ( kind = 4 ) va
  integer ( kind = 4 ) vb
  integer ( kind = 4 ) vc
  integer ( kind = 4 ) vm(*)
  real    ( kind = 8 ) vcl(3,*)
!
!  TOPN is index of top of stack of non-visible boundary faces.
!  TOPT is index of top of stack of boundary faces to be tested.
!
  topn = topnv
  topt = 0
  fc(7,topv) = 0
  k = -fc(5,topv)

  do j = 1, 3
    nbr = bf(j,k)
    if ( fc(7,nbr) == -1 ) then
      fc(7,nbr) = topt
      topt = nbr
    end if
  end do

  do while ( topt /= 0 )

    ptr = topt
    topt = fc(7,ptr)
    va = vm(fc(1,ptr))
    vb = vm(fc(2,ptr))
    vc = vm(fc(3,ptr))
    op = opside ( vcl(1,va), vcl(1,vb), vcl(1,vc), ctr, pt )

    if ( op == 2 ) then
      ierr = 301
      return
    end if

    if ( op == 1 ) then

      fc(7,ptr) = topv
      topv = ptr
      k = -fc(5,ptr)

      do j = 1, 3
        nbr = bf(j,k)
        if ( fc(7,nbr) == -1 ) then
          fc(7,nbr) = topt
          topt = nbr
        end if
      end do

    else

      fc(7,ptr) = topn
      topn = ptr

    end if

  end do
!
!  For boundary faces not visible from PT, set FC(7,*) = -1.
!
  do while ( topn /= 0 )
    ptr = topn
    topn = fc(7,ptr)
    fc(7,ptr) = -1
  end do

  return
end
