subroutine vbedg ( x, y, vcl, til, tnbr, ltri, ledg, rtri, redg )

!*****************************************************************************80
!
!! VBEDG determines the boundary edges of a 2D triangulation.
!
!  Discussion:
!
!    This routine determines boundary edges of a 2D triangulation which are
!    visible from point (X,Y) outside convex hull.
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
!    Input, real ( kind = 8 ) X, Y, the 2D point outside convex hull.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) TIL(1:3,1:*), the triangle incidence list.
!
!    Input, integer ( kind = 4 ) TNBR(1:3,1:*), the triangle neighbor list; negative
!    values are used for links of counterclockwise linked list of boundary
!    edges; LINK = -(3*I + J-1) where I, J = triangle, edge index.
!
!    Input, integer ( kind = 4 ) LTRI, LEDG, if LTRI /= 0 then they are assumed to be as
!    defined below and are not changed, else they are updated.  LTRI is
!    the index of boundary triangle to left of leftmost boundary
!    triangle visible from (X,Y).  LEDG is the boundary edge of triangle
!    LTRI to left of leftmost boundary edge visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of boundary triangle
!    to begin search at.  On output, the index of rightmost boundary triangle
!    visible from (X,Y)
!
!    Input/output, integer ( kind = 4 ) REDG.  On input, the edge of triangle RTRI that is
!    visible from (X,Y).  On output, the edge of triangle RTRI that is
!    visible from (X,Y)
!
  implicit none

  integer ( kind = 4 ) a
  integer ( kind = 4 ) b
  integer ( kind = 4 ) e
  integer ( kind = 4 ) l
  logical              ldone
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,*)
  integer ( kind = 4 ) tnbr(3,*)
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
!
!  Find rightmost visible boundary edge using links, then possibly
!  leftmost visible boundary edge using triangle neighbor info.
!
  if ( ltri == 0 ) then
    ldone = .false.
    ltri = rtri
    ledg = redg
  else
    ldone = .true.
  end if

  do

    l = -tnbr(redg,rtri)
    t = l / 3
    e = mod ( l, 3 ) + 1
    a = til(e,t)

    if ( e <= 2 ) then
      b = til(e+1,t)
    else
      b = til(1,t)
    end if

    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0D+00)

    if ( lr <= 0 ) then
      exit
    end if

    rtri = t
    redg = e

  end do

  if ( ldone ) then
    return
  end if

  t = ltri
  e = ledg

  do

    b = til(e,t)

    if ( 2 <= e ) then
      e = e - 1
    else
      e = 3
    end if

    do while ( 0 < tnbr(e,t) )

      t = tnbr(e,t)

      if ( til(1,t) == b ) then
        e = 3
      else if ( til(2,t) == b ) then
        e = 1
      else
        e = 2
      end if

    end do

    a = til(e,t)
    lr = lrline(x,y,vcl(1,a),vcl(2,a),vcl(1,b),vcl(2,b),0.0D+00)

    if ( lr <= 0 ) then
      exit
    end if

  end do

  ltri = t
  ledg = e

  return
end
