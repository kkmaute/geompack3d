subroutine fndsep ( angac1, xr, yr, nvrt, xc, yc, ivis, theta, nv, iv,  &
  vcl, pvl, iang, angsep, i1, i2, wkang )

!*****************************************************************************80
!
!! FNDSEP finds separators to resolve a reflex vertex.
!
!  Discussion:
!
!    This routine finds 1 or 2 separators which can resolve a reflex vertex
!    (XR,YR) using a max-min angle criterion from list of vertices
!    in increasing polar angle with respect to the reflex vertex.
!
!    Preference is given to 1 separator.
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
!    Input, real ( kind = 8 ) ANGAC1, the angle tolerance parameter used
!    for preference in accepting one separator.
!
!    Input, real ( kind = 8 ) XR, YR, the coordinates of reflex vertex.
!
!    Input, integer ( kind = 4 ) NVRT, (number of vertices) - 1.
!
!    Input, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT), the vertex coordinates
!    of possible endpoints of a separator.
!
!    Input, integer ( kind = 4 ) IVIS(0:NVRT), contains information about the vertices of
!    XC, YC arrays with respect to the polygon vertex list; if
!    IVIS(I) greater than 0 then vertex (XC(I),YC(I)) has index IVIS(I)
!    in PVL; if IVIS(I) < 0 then vertex (XC(I),YC(I)) is on
!    the edge joining vertices with indices -IVIS(I) and
!    SUCC(-IVIS(I)) in PVL.
!
!    Input, real ( kind = 8 ) THETA(0:NVRT), the polar angles of vertices
!    in increasing order; THETA(NVRT) is the interior angle of reflex vertex;
!    THETA(I), 0 <= I, is the polar angle of (XC(I),YC(I))
!    with respect to reflex vertex.
!
!    Input, integer ( kind = 4 ) NV, (number of vertices to be considered as endpoint of a
!    separator) - 1.
!
!    Input, integer ( kind = 4 ) IV(0:NV), the indices of vertices in XC, YC arrays to be
!    considered as endpoint of a separator; angle between
!    consecutive vertices is assumed to be < 180 degrees.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) PVL(1:4,1:*), real ( kind = 8 ) IANG(1:*), the polygon
!    vertex list, interior angles.
!
!    Output, real ( kind = 8 ) ANGSEP, the minimum of the 4 or 7 angles at the
!    boundary resulting from 1 or 2 separators, respectively.
!
!    Output, integer ( kind = 4 ) I1, I2, the indices of endpoints of separators in XC,
!    YC arrays; I2 = -1 if there is only one separator, else I1 < I2.
!
!    Workspace, real ( kind = 8 ) WKANG(0:NV).
!
  implicit none

  real    ( kind = 8 ) ang
  real    ( kind = 8 ) angac1
  real    ( kind = 8 ) angsep
  real    ( kind = 8 ) angsp2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  real    ( kind = 8 ) iang(*)
  integer ( kind = 4 ) ii
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvrt
  integer ( kind = 4 ) iv(0:nv)
  integer ( kind = 4 ) ivis(0:nvrt)
  real    ( kind = 8 ) minang
  integer ( kind = 4 ) p
  real    ( kind = 8 ) phi
  real    ( kind = 8 ), parameter :: pi = 3.141592653589793D+00
  integer ( kind = 4 ) pvl(4,*)
  integer ( kind = 4 ) q
  integer ( kind = 4 ) r
  real    ( kind = 8 ) theta(0:nvrt)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(2,*)
  real    ( kind = 8 ) wkang(0:nv)
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xr
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) yr

  tol = 100.0D+00 * epsilon ( tol )
!
!  Determine the vertices in the inner cone - indices P to Q.
!
  i = 0
  p = -1
  phi = theta(nvrt) - pi + tol

  do while ( p < 0 )

    if ( phi <= theta(iv(i)) ) then
      p = i
    else
      i = i + 1
    end if

  end do

  i = nv
  q = -1
  phi = pi - tol

  do while ( q < 0 )

    if ( theta(iv(i)) <= phi ) then
      q = i
    else
      i = i - 1
    end if

  end do
!
!  Use the max-min angle criterion to find the best separator
!  in inner cone.
!
  angsep = 0.0D+00

  do i = p, q

    k = iv(i)
    ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
      vcl, pvl, iang )

    if ( angsep < ang ) then
      angsep = ang
      ii = iv(i)
    end if

  end do

  angsp2 = angsep
  if ( angac1 <= angsep ) then
    go to 110
  end if
!
!  If the best separator in inner cone is not 'good' enough,
!  use max-min angle criterion to try to find a better pair
!  of separators from the right and left cones.
!
  nr = 0
  nl = 0

  do r = 0, p-1

    wkang(r) = 0.0D+00

    if ( angsep < theta(iv(r)) ) then

      k = iv(r)

      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( angsep < ang ) then
        nr = nr + 1
        wkang(r) = ang
      end if

    end if

  end do

  if ( nr == 0 ) then
    go to 110
  end if

  phi = theta(nvrt) - angsep

  do l = q+1, nv

    wkang(l) = 0.0D+00

    if ( theta(iv(l)) < phi ) then

      k = iv(l)
      ang = minang ( xr, yr, xc(k), yc(k), ivis(k), theta(k), theta(nvrt), &
        vcl, pvl, iang )

      if ( angsep < ang ) then
        nl = nl + 1
        wkang(l) = ang
      end if

    end if

  end do

  if ( nl == 0 ) then
    go to 110
  end if
!
!  Check all possible pairs for the best pair of separators
!  in the right and left cones.
!
  m = nv

  do r = p-1, 0, -1

     if ( q < m .and. angsp2 < wkang(r) ) then

      phi = theta(iv(r))

80    continue

      if ( q < m .and. &
        ( wkang(m) <= angsp2 .or. pi - tol < theta(iv(m)) - phi ) ) then
        m = m - 1
        go to 80
      end if

      do l = q+1, m

         if ( angsp2 < wkang(l) ) then

            ang = min ( theta(iv(l)) - phi, wkang(r), wkang(l) )

            if ( angsp2 < ang ) then
             angsp2 = ang
             i1 = iv(r)
             i2 = iv(l)
            end if

         end if

       end do

     end if

  end do
!
!  Choose 1 or 2 separators based on max-min angle criterion or
!  ANGAC1 parameter.
!
  110 continue

  if ( angsp2 <= angsep ) then
    i1 = ii
    i2 = -1
  else
    angsep = angsp2
  end if

  return
end
