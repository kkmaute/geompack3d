subroutine intpg ( nvrt, xc, yc, ctrx, ctry, arpoly, hflag, umdf, wsq, nev, &
  ifv, listev, ivrt, edgval, vrtval, vcl, mdfint, mean, stdv, mdftr )

!*****************************************************************************80
!
!! INTPG integrates a mesh distribution function over a polygon.
!
!  Discussion:
!
!    This routine computes the integral of MDF2(X,Y) [heuristic mdf] or
!    UMDF(X,Y) [user-supplied mdf] in a convex polygon.
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
!    Input, integer ( kind = 4 ) NVRT, the number of vertices in polygon.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the coordinates of
!    polygon vertices in counterclockwise order, translated so that centroid is
!    at origin; (XC(0),YC(0)) = (XC(NVRT),YC(NVRT)).
!
!    Input, real ( kind = 8 ) CTRX, CTRY, the coordinates of the centroid
!    before translation.
!
!    Input, real ( kind = 8 ) ARPOLY, the area of polygon.
!
!    Input, logical HFLAG, is TRUE if heuristic mdf, FALSE if
!    user-supplied mdf.
!
!    Input, external real ( kind = 8 ) UMDF(X,Y), the user-supplied mdf
!    with d.p arguments.
!
!    Input, real ( kind = 8 ) WSQ, square of width of original polygon
!    of decomposition.
!
!    Input, integer ( kind = 4 ) NEV, IFV, LISTEV(1:NEV), the output from routine PRMDF2.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), real ( kind = 8 ) EDGVAL(1:*), VRTVAL(1:*),
!    arrays output from DSMDF2; if .NOT. HFLAG then only first array exists.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Output, real ( kind = 8 ) MDFINT, the integral of mdf in polygon.
!
!    Output, real ( kind = 8 ) MEAN, the mean mdf value in polygon.
!
!    Output, real ( kind = 8 ) STDV, the standard deviation of mdf in polygon.
!
!    Output, real ( kind = 8 ) MDFTR(0:NVRT-1), mean mdf value in each
!    triangle of polygon; triangles are determined by polygon vertices
!    and centroid.
!
  implicit none

  integer ( kind = 4 ) nev
  integer ( kind = 4 ), parameter :: nqpt = 3
  integer ( kind = 4 ) nvrt

  real    ( kind = 8 ) areatr
  real    ( kind = 8 ) arpoly
  real    ( kind = 8 ) ctrx
  real    ( kind = 8 ) ctry
  real    ( kind = 8 ) d
  real    ( kind = 8 ) edgval(*)
  logical              hflag
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev(nev)
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mdfint
  real    ( kind = 8 ) mdfsqi
  real    ( kind = 8 ) mdftr(0:nvrt-1)
  real    ( kind = 8 ) mean
  real    ( kind = 8 ), parameter, dimension (3,nqpt) :: qc = reshape ( &
    (/ 0.6666666666666666D+00, 0.1666666666666667D+00, &
       0.1666666666666667D+00, 0.1666666666666667D+00, &
       0.6666666666666666D+00, 0.1666666666666667D+00, &
       0.1666666666666667D+00, 0.1666666666666667D+00, &
       0.6666666666666666D+00/), (/ 3, nqpt /) )
  real    ( kind = 8 ) s
  real    ( kind = 8 ) stdv
  real    ( kind = 8 ) sum1
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) umdf
  real    ( kind = 8 ) val
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) wsq
  real    ( kind = 8 ), parameter, dimension ( nqpt ) :: wt = (/ &
    0.3333333333333333D+00, 0.3333333333333333D+00, 0.3333333333333333D+00 /)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) x0
  real    ( kind = 8 ) x1
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xx
  real    ( kind = 8 ) y
  real    ( kind = 8 ) y0
  real    ( kind = 8 ) y1
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) yy
!
!  NQPT is number of quad points for numerical integration in triangle
!  WT(I) is weight of Ith quadrature point
!  QC(1:3,I) are barycentric coordinates of Ith quadrature point
!
  mdfint = 0.0D+00
  mdfsqi = 0.0D+00

  do l = 0, nvrt-1

    areatr = 0.5D+00 * ( xc(l) * yc(l+1) - xc(l+1) * yc(l) )
    sum1 = 0.0D+00
    sum2 = 0.0D+00

    do m = 1, nqpt

      xx = qc(1,m) * xc(l) + qc(2,m) * xc(l+1)
      yy = qc(1,m) * yc(l) + qc(2,m) * yc(l+1)
!
!  Insert code for function MDF2 to reduce number of calls.
!
      if ( hflag ) then

        x = xx + ctrx
        y = yy + ctry
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

              if ( 0.0D+00 <= x0 * x1 + y0 * y1 ) then
                d = x0**2 + y0**2
              else
                d = (x1*y0 - y1*x0)**2 / (x1**2 + y1**2)
              end if

            end if

            d = max ( 0.25D+00 * d, edgval(k) )
            s = min ( s, d )

          end if

        end do

        val = 1.0D+00/s

      else

        val = umdf ( xx+ctrx, yy+ctry )

      end if

      temp = wt(m) * val
      sum1 = sum1 + temp
      sum2 = sum2 + temp*val

    end do

    mdftr(l) = sum1
    mdfint = mdfint + sum1 * areatr
    mdfsqi = mdfsqi + sum2 * areatr

  end do

  mean = mdfint / arpoly
  stdv = mdfsqi / arpoly - mean**2
  stdv = sqrt ( max ( stdv, 0.0D+00 ) )

  return
end
