subroutine swapec ( i, top, maxst, btri, bedg, vcl, til, tnbr, stack, ierr )

!*****************************************************************************80
!
!! SWAPEC swaps diagonal edges in a 2D triangulation
!
!  Discussion:
!
!    This routine swaps diagonal edges in a 2D triangulation based on empty
!    circumcircle criterion until all triangles are Delaunay, given
!    that I is index of new vertex added to triangulation.
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
!    Input, integer ( kind = 4 ) I, the index in VCL of new vertex.
!
!    Input/output, integer ( kind = 4 ) TOP, the index of top of stack, must be greater
!    than or equal to 0.
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array.
!
!    Input/output, integer ( kind = 4 ) BTRI, BEDG, if positive, these are triangle and
!    edge index of a boundary edge whose updated indices must be recorded.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input/output, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input/output, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) STACK(1:TOP), the index of initial triangles (involving
!    vertex I) put in stack; the edges opposite I should be in interior.
!
!    Workspace, integer STACK(TOP+1:MAXST), the used as stack.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxst

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) bedg
  integer ( kind = 4 ) btri
  integer ( kind = 4 ) c
  integer ( kind = 4 ) diaedg
  integer ( kind = 4 ) e
  integer ( kind = 4 ) ee
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) f
  integer ( kind = 4 ) fm1
  integer ( kind = 4 ) fp1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) l
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) r
  integer ( kind = 4 ) s
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) swap
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  integer ( kind = 4 ) top
  integer ( kind = 4 ) tt
  integer ( kind = 4 ) u
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) vcl(2,*)
!
!  Determine whether triangles in stack are Delaunay, and swap
!  diagonal edge of convex quadrilateral if not.
!
  ierr = 0
  x = vcl(1,i)
  y = vcl(2,i)

10 continue

  if ( top <= 0 ) then
    return
  end if

  t = stack(top)
  top = top - 1

  if ( til(1,t) == i ) then
    e = 2
    b = til(3,t)
  else if ( til(2,t) == i ) then
    e = 3
    b = til(1,t)
  else
    e = 1
    b = til(2,t)
  end if

  a = til(e,t)
  u = tnbr(e,t)

  if ( tnbr(1,u) == t ) then
    f = 1
    c = til(3,u)
  else if ( tnbr(2,u) == t ) then
    f = 2
    c = til(1,u)
  else
    f = 3
    c = til(2,u)
  end if

  swap = diaedg(x,y,vcl(1,a),vcl(2,a),vcl(1,c),vcl(2,c),vcl(1,b), &
    vcl(2,b))

  if ( swap == 1 ) then

    em1 = e - 1
    if ( em1 == 0) em1 = 3
    ep1 = e + 1
    if ( ep1 == 4) ep1 = 1
    fm1 = f - 1
    if ( fm1 == 0) fm1 = 3
    fp1 = f + 1
    if ( fp1 == 4) fp1 = 1
    til(ep1,t) = c
    til(fp1,u) = i
    r = tnbr(ep1,t)
    s = tnbr(fp1,u)
    tnbr(ep1,t) = u
    tnbr(fp1,u) = t
    tnbr(e,t) = s
    tnbr(f,u) = r

    if ( 0 < tnbr(fm1,u) ) then
      top = top + 1
      stack(top) = u
    end if

    if ( 0 < s ) then

      if ( tnbr(1,s) == u ) then
        tnbr(1,s) = t
      else if ( tnbr(2,s) == u ) then
        tnbr(2,s) = t
      else
        tnbr(3,s) = t
      end if

      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        return
      end if

      stack(top) = t

    else

      if ( u == btri .and. fp1 == bedg ) then
        btri = t
        bedg = e
      end if

      l = -( 3 * t + e - 1 )
      tt = t
      ee = em1

20    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == a ) then
          ee = 3
        else if ( til(2,tt) == a ) then
          ee = 1
        else
          ee = 2
        end if

        go to 20

      end if

      tnbr(ee,tt) = l

    end if

    if ( 0 < r ) then

      if ( tnbr(1,r) == t ) then
        tnbr(1,r) = u
      else if ( tnbr(2,r) == t ) then
        tnbr(2,r) = u
      else
        tnbr(3,r) = u
      end if

    else

      if ( t == btri .and. ep1 == bedg ) then
        btri = u
        bedg = f
      end if

      l = -(3*u + f-1)
      tt = u
      ee = fm1

30    continue

      if ( 0 < tnbr(ee,tt) ) then

        tt = tnbr(ee,tt)

        if ( til(1,tt) == b ) then
          ee = 3
        else if ( til(2,tt) == b ) then
          ee = 1
        else
          ee = 2
        end if

        go to 30

      end if

      tnbr(ee,tt) = l

    end if

    if ( msglvl == 4 ) then
      write ( *,600) 2,vcl(1,a),vcl(2,a), &
           vcl(1,b),vcl(2,b),x,y,vcl(1,c),vcl(2,c)
    end if

  end if

  go to 10

  600 format (1x,i7,4f15.7/8x,4f15.7)

end
