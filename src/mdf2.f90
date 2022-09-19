function mdf2 ( x, y, wsq, nev, ifv, listev, ivrt, edgval, vrtval, vcl )

!*****************************************************************************80
!
!! MDF2 evaluates a heuristic mesh distribution function in 2D.
!
!  Discussion:
!
!    This routine evaluates a heuristic mesh distribution function at (X,Y).
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
!    Input, real ( kind = 8 ) X, Y, the coordinates of point.
!
!    Input, real ( kind = 8 ) WSQ, the square of width of polygon
!    containing (X,Y).
!
!    Input, integer ( kind = 4 ) NEV, IFV, LISTEV(1:NEV), output from routine PRMDF2.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), real ( kind = 8 ) EDGVAL(1:*),
!    real ( kind = 8 ) VRTVAL(1:*), arrays output from DSMDF2.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list
!
!    Output, real ( kind = 8 ) MDF2, the reciprocal of square of length
!    scale at (X,Y).
!
  implicit none

  integer ( kind = 4 ) nev

  real    ( kind = 8 ) d
  real    ( kind = 8 ) edgval(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) listev(nev)
  real    ( kind = 8 ) mdf2
  real    ( kind = 8 ) s
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) wsq
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) y1

  s = wsq

  do i = 1, nev

    k = listev(i)

    if ( k < 0 ) then

      k = -k
      d = (vcl(1,k) - x)**2 + (vcl(2,k) - y)**2
      d = max ( 0.25D+00 * d, vrtval(k) )
      s = min(s,d)

    else

      kp1 = k + 1
      if ( i == nev .and. 0 < ifv ) then
        kp1 = ifv
      end if

      j = ivrt(kp1)
      x0 = x - vcl(1,j)
      y0 = y - vcl(2,j)
      x1 = vcl(1,ivrt(k)) - vcl(1,j)
      y1 = vcl(2,ivrt(k)) - vcl(2,j)

      if ( x0*x1 + y0*y1 <= 0.0D+00 ) then
        d = x0**2 + y0**2
      else
        x0 = x0 - x1
        y0 = y0 - y1
        if ( 0.0D+00 <= x0*x1 + y0*y1 ) then
          d = x0**2 + y0**2
        else
          d = (x1*y0 - y1*x0)**2/(x1**2 + y1**2)
        end if
      end if

      d = max ( 0.25D+00 * d, edgval(k) )
      s = min(s,d)

    end if

  end do

  mdf2 = 1.0D+00/s

  return
end
