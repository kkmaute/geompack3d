subroutine prmdf2 ( ipoly, wsq, ivrt, xivrt, edgval, vrtval, nev, ifv, &
  listev )

!*****************************************************************************80
!
!! PRMDF2 does preprocessing for the mesh distribution function evaluation.
!
!  Discussion:
!
!    This routine is a preprocessing step for evaluating a mesh distribution
!    function in polygon IPOLY - the edges and vertices for
!    which distances must be computed are determined.
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
!    Input, integer ( kind = 4 ) IPOLY, the index of polygon.
!
!    Input, real ( kind = 8 ) WSQ, the square of width of polygon IPOLY.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), the indices of polygon vertices in VCL,
!    ordered by polygon.
!
!    Input, integer ( kind = 4 ) XIVRT(1:*), pointer to first vertex of each polygon in
!    IVRT; vertices of polygon IPOLY are IVRT(I) for I from XIVRT(IPOLY)
!    to XIVRT(IPOLY+1)-1.
!
!    Input, real ( kind = 8 ) EDGVAL(1:*), the value associated with each
!    edge of decomposition.
!
!    Input, real ( kind = 8 ) VRTVAL(1:*), the value associated with each
!    vertex of decomposition.
!
!    Output, integer ( kind = 4 ) NEV, the number of edges and vertices for which distances
!    must be evaluated.
!
!    Output, integer ( kind = 4 ) IFV, the index of first vertex XIVRT(IPOLY) if LISTEV(NEV)
!    = XIVRT(IPOLY+1) - 1; 0 otherwise.
!
!    Output, integer ( kind = 4 ) LISTEV(1:*), array of length
!    <= [XIVRT(IPOLY+1)-XIVRT(IPOLY)]*2 containing indices of edges and
!    vertices mentioned above; indices of vertices are negated.
!
  implicit none

  real    ( kind = 8 ) edgval(*)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ifv
  integer ( kind = 4 ) im1
  integer ( kind = 4 ) ipoly
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ) nev
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) wsq
  integer ( kind = 4 ) xivrt(*)

  ifv = 0
  nev = 0
  im1 = xivrt(ipoly+1) - 1
  l = im1

  do i = xivrt(ipoly), l

    j = ivrt(i)

    if ( vrtval(j) < min ( edgval(i), edgval(im1) ) ) then
      nev = nev + 1
      listev(nev) = -j
    end if

    if ( edgval(i) < wsq ) then
      nev = nev + 1
      listev(nev) = i
    end if

    im1 = i

  end do

  if ( 0 < nev ) then
    if ( listev(nev) == l ) then
      ifv = xivrt(ipoly)
    end if
  end if

  return
end
