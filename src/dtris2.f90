subroutine dtris2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )

!*****************************************************************************80
!
!! DTRIS2 constructs the Delaunay triangulation of vertices in 2D.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of 2D vertices using
!    incremental approach and diagonal edge swaps.  Vertices are
!    first sorted in lexicographically increasing (x,y) order, and
!    then are inserted one at a time from outside the convex hull.
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
!    Input, integer ( kind = 4 ) NPT, the number of 2D points (vertices).
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array;
!    should be about NPT to be safe, but MAX ( 10, 2*LOG2(NPT) ) usually enough.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), coordinates of 2D vertices.
!
!    Input/output, IND(1:NPT), the indices in VCL of vertices to be
!    triangulated.  On output, permuted due to sorting.
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; elements
!    are indices of VCL; vertices of triangles are in counterclockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list; positive
!    elements are indices of TIL; negative elements are used for links
!    of counterclockwise linked list of boundary edges;
!    LINK = -(3*I + J-1) where I, J = triangle, edge index;
!    TNBR(J,I) refers to the neighbor along edge from vertex
!    J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), used for stack of triangles for which
!    circumcircle test must be made.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxst
  integer ( kind = 4 ) npt

  real    ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ind(npt)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) l
  integer ( kind = 4 ) ledg
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) ltri
  integer ( kind = 4 ) m
  integer ( kind = 4 ) m1
  integer ( kind = 4 ) m2
  integer ( kind = 4 ), parameter :: msglvl = 0
  integer ( kind = 4 ) n
  integer ( kind = 4 ) ntri
  integer ( kind = 4 ) redg
  integer ( kind = 4 ) rtri
  integer ( kind = 4 ) stack(maxst)
  integer ( kind = 4 ) t
  integer ( kind = 4 ) til(3,npt*2)
  integer ( kind = 4 ) tnbr(3,npt*2)
  real    ( kind = 8 ) tol
  integer ( kind = 4 ) top
  real    ( kind = 8 ) vcl(2,*)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
!
!  Sort vertices by increasing (x,y).
!
  call dhpsrt ( 2, npt, 2, vcl, ind )
!
!  Check that the data is not degenerate.
!
  m1 = ind(1)

  do i = 2, npt

    m = m1
    m1 = ind(i)

    do j = 1, 2
      cmax = max ( abs(vcl(j,m)), abs(vcl(j,m1)) )
      if ( tol * cmax < abs(vcl(j,m) - vcl(j,m1)) .and. tol < cmax ) then
        go to 20
      end if
    end do

    ierr = 224
    return

20  continue

  end do
!
!  Staring with points M1 and M2, find the first point M that is
!  "reasonably" non-colinear.
!
  m1 = ind(1)
  m2 = ind(2)
  j = 3

  do

    if ( npt < j ) then
      ierr = 225
      return
    end if

    m = ind(j)
    lr = lrline ( vcl(1,m), vcl(2,m), vcl(1,m1), vcl(2,m1), vcl(1,m2), &
      vcl(2,m2), 0.0D+00 )

    if ( lr /= 0 ) then
      exit
    end if

    j = j + 1

  end do
!
!  Set up the initial triangle information for M1, M2 and M (and any
!  in-between points we may have skipped over while searching for M.
!
  ntri = j - 2

  if ( lr == -1 ) then

    til(1,1) = m1
    til(2,1) = m2
    til(3,1) = m
    tnbr(3,1) = -3

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m1
      til(2,i) = m2
      til(3,i) = m
      tnbr(1,i-1) = -3*i
      tnbr(2,i-1) = i
      tnbr(3,i) = i - 1
    end do

    tnbr(1,ntri) = -3*ntri - 1
    tnbr(2,ntri) = -5
    ledg = 2
    ltri = ntri

  else

    til(1,1) = m2
    til(2,1) = m1
    til(3,1) = m
    tnbr(1,1) = -4

    do i = 2, ntri
      m1 = m2
      m2 = ind(i+1)
      til(1,i) = m2
      til(2,i) = m1
      til(3,i) = m
      tnbr(3,i-1) = i
      tnbr(1,i) = -3*i - 3
      tnbr(2,i) = i - 1
    end do

    tnbr(3,ntri) = -3*ntri
    tnbr(2,1) = -3*ntri - 2
    ledg = 2
    ltri = 1

  end if

  if ( msglvl == 4 ) then
    m2 = ind(1)
    write ( *,'(i7,4f15.7)') 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    do i = 2, j-1
      m1 = m2
      m2 = ind(i)
      write ( *,'(i7,4f15.7)') 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
      write ( *,'(i7,4f15.7)') 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    end do
  end if
!
!  Insert vertices one at a time from outside convex hull, determine
!  visible boundary edges, and apply diagonal edge swaps until
!  Delaunay triangulation of vertices (so far) is obtained.
!
  top = 0

  do i = j+1, npt

    if ( msglvl == 4 ) then
      write ( *,'(i7,4f15.7)') i
    end if

    m = ind(i)
    m1 = til(ledg,ltri)

    if ( ledg <= 2 ) then
      m2 = til(ledg+1,ltri)
    else
      m2 = til(1,ltri)
    end if

    lr = lrline(vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1),vcl(1,m2), &
      vcl(2,m2),0.0D+00)

    if ( 0 < lr ) then
      rtri = ltri
      redg = ledg
      ltri = 0
    else
      l = -tnbr(ledg,ltri)
      rtri = l / 3
      redg = mod(l,3) + 1
    end if

    call vbedg(vcl(1,m),vcl(2,m),vcl,til,tnbr,ltri,ledg,rtri,redg)
    n = ntri + 1
    l = -tnbr(ledg,ltri)

    do

      t = l / 3
      e = mod(l,3) + 1
      l = -tnbr(e,t)
      m2 = til(e,t)

      if ( e <= 2 ) then
        m1 = til(e+1,t)
      else
        m1 = til(1,t)
      end if

      ntri = ntri + 1
      tnbr(e,t) = ntri
      til(1,ntri) = m1
      til(2,ntri) = m2
      til(3,ntri) = m
      tnbr(1,ntri) = t
      tnbr(2,ntri) = ntri - 1
      tnbr(3,ntri) = ntri + 1
      top = top + 1

      if ( maxst < top ) then
        ierr = 8
        return
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,'(i7,4f15.7)') 1,vcl(1,m),vcl(2,m), vcl(1,m2),vcl(2,m2)
      end if

      if ( t == rtri .and. e == redg ) then
        exit
      end if

    end do

    if ( msglvl == 4 ) then
      write ( *,'(i7,4f15.7)') 1,vcl(1,m),vcl(2,m), vcl(1,m1),vcl(2,m1)
    end if

    tnbr(ledg,ltri) = -3*n - 1
    tnbr(2,n) = -3*ntri - 2
    tnbr(3,ntri) = -l
    ltri = n
    ledg = 2
    call swapec(m,top,maxst,ltri,ledg,vcl,til,tnbr,stack, ierr )

    if ( ierr /= 0 ) then
      return
    end if

  end do

  if ( msglvl == 4 ) then
    write ( *, '(i7,4f15.7)' ) npt+1
  end if

  return
end
