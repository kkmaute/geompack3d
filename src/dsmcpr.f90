subroutine dsmcpr ( nhole, nvbc, vcl, maxhv, maxpv, maxho, nvc, npolg, &
  nvert, nhola, regnum, hvl, pvl, iang, holv, ierr )

!*****************************************************************************80
!
!! DSMCPR initializes the polygonal decomposition data structure.
!
!  Discussion:
!
!    This routine initializes the polygonal decomposition data structure
!    given a multiply-connected polygonal region with 1 outer
!    boundary curve and 0 or more inner boundary curves of holes.
!
!  Author:
!
!    Barry Joe,
!    Department of Computing Science,
!    University of Alberta,
!    Edmonton, Alberta, Canada  T6G 2H1
!    Email: barry@cs.ualberta.ca
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
!    Input, integer ( kind = 4 ) NHOLE, the number of holes in region.
!
!    Input, integer ( kind = 4 ) NVBC(1:NHOLE+1), the number of vertices per boundary
!    curve; first boundary curve is the outer boundary of the region.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:NVC), vertex coordinates of boundary
!    curves in counterclockwise order; NVC = NVBC(1) + ... + NVBC(NHOLE+1);
!    positions 1 to NVBC(1) of VCL contain the vertex coordinates of the
!    outer boundary in counterclockwise order; positions NVBC(1)+1 to
!    NVBC(1)+NVBC(2) contain the vertex coordinates of the
!    first hole boundary in counterclockwise order, etc.
!
!    Input, integer ( kind = 4 ) MAXHV, the maximum size available for HVL, REGNUM arrays,
!    should be greater than or equal to NHOLE + 1.
!
!    Input, integer ( kind = 4 ) MAXPV, the maximum size available for PVL, IANG arrays;
!    should be greater than or equal to NVC.
!
!    Input, integer ( kind = 4 ) MAXHO, the maximum size available for HOLV array; should be
!    greater than or equal to NHOLE*2.
!
!    Output, integer ( kind = 4 ) NVC, the number of vertex coordinates, set to
!    sum of NVBC(I).
!
!    Output, integer ( kind = 4 ) NPOLG, the number of polygonal subregions, set to 1.
!
!    Output, integer ( kind = 4 ) NVERT, the number of vertices in PVL, set to NVC.
!
!    Output, integer ( kind = 4 ) NHOLA, the number of attached holes, set to 0.
!
!    Output, integer ( kind = 4 ) REGNUM(1:1), the region number of only subregion, set to 1.
!
!    [Note: Above 4 parameters are for consistency with DSPGDC.]
!
!    Output, integer ( kind = 4 ) HVL(1:NHOLE+1), the head vertex list; first entry
!    is the head vertex (index in PVL) of outer boundary curve; next
!    NHOLE entries contain the head vertex of a hole.
!
!    Output, integer ( kind = 4 ) PVL(1:4,1:NVC),IANG(1:NVC), the polygon vertex list
!    and interior angles; vertices of outer boundary curve are in
!    counterclockwise order followed by vertices of each hole in CW hole;
!    vertices of each polygon are in a circular linked list; see
!    routine DSPGDC for more details of this data structure.
!
!    Output, integer ( kind = 4 ) HOLV(1:NHOLE*2), the indices in PVL of top and bottom
!    vertices of holes; first (last) NHOLE entries are for top (bottom)
!    vertices; top (bottom) vertices are sorted in decreasing
!    (increasing) lexicographic (y,x) order of coordinates.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxho
  integer ( kind = 4 ) maxpv
  integer ( kind = 4 ) nhole

  real    ( kind = 8 ) angle
  integer ( kind = 4 ), parameter :: edgv = 4
  integer ( kind = 4 ) holv(maxho)
  integer ( kind = 4 ) hvl(nhole+1)
  integer ( kind = 4 ) i
  real    ( kind = 8 ) iang(maxpv)
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) iv
  integer ( kind = 4 ) ivs
  integer ( kind = 4 ) j
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) lv
  integer ( kind = 4 ) lvp
  integer ( kind = 4 ) lvs
  integer ( kind = 4 ) maxhv
  integer ( kind = 4 ) nhola
  integer ( kind = 4 ) npolg
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvbc(nhole+1)
  integer ( kind = 4 ) nvc
  integer ( kind = 4 ) nvert
  integer ( kind = 4 ) nvs
  integer ( kind = 4 ), parameter :: polg = 2
  integer ( kind = 4 ) pvl(4,maxpv)
  integer ( kind = 4 ) regnum(1)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) vcl(2,*)

  ierr = 0
  nvc = sum ( nvbc(1:nhole+1) )

  npolg = 1
  nvert = nvc
  nhola = 0
  regnum(1) = 1

  if ( maxhv < nhole + 1 ) then
    ierr = 4
    return
  else if ( maxpv < nvc ) then
    ierr = 5
    return
  else if ( maxho < nhole + nhole ) then
    ierr = 2
    return
  end if
!
!  Initialize HVL, PVL arrays.
!
20 continue

  hvl(1) = 1
  nv = nvbc(1)

  do i = 1, nv
    pvl(loc,i) = i
    pvl(polg,i) = 1
    pvl(succ,i) = i + 1
    pvl(edgv,i) = 0
  end do

  pvl(succ,nv) = 1

  do j = 1, nhole
    hvl(j+1) = nv + 1
    nvs = nv + nvbc(j+1)
    do i = nv+1, nvs
      pvl(loc,i) = i
      pvl(polg,i) = 1
      pvl(succ,i) = i - 1
      pvl(edgv,i) = 0
    end do
    pvl(succ,nv+1) = nvs
    nv = nvs
  end do
!
!  Initialize IANG array.
!
  do i = 1, nhole+1

    j = hvl(i)
    lvp = pvl(loc,j)
    iv = pvl(succ,j)
    lv = pvl(loc,iv)

    do

      ivs = pvl(succ,iv)
      lvs = pvl(loc,ivs)
      iang(iv) = angle ( vcl(1,lvp), vcl(2,lvp), vcl(1,lv), vcl(2,lv), &
        vcl(1,lvs), vcl(2,lvs) )

      if ( iv == j ) then
        exit
      end if

      lvp = lv
      iv = ivs
      lv = lvs

    end do

  end do
!
!  Initialize HOLV array.
!
  if ( 0 < nhole ) then
    call holvrt ( nhole, vcl, hvl(2), pvl, holv )
  end if

  return
end
