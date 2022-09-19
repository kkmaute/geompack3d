subroutine dtriw2 ( npt, maxst, vcl, ind, ntri, til, tnbr, stack, ierr )

!*****************************************************************************80
!
!! DTRIW2 constructs a Delaunay triangulation of vertices in 2D.
!
!  Discussion:
!
!    This routine constructs a Delaunay triangulation of 2D vertices using
!    incremental approach and diagonal edge swaps. Vertices are
!    inserted one at a time in order given by IND array. The initial
!    triangles created due to a new vertex are obtained by a walk
!    through the triangulation until location of vertex is known.
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
!    Input, integer ( kind = 4 ) NPT, the number of 2D points (vertices).
!
!    Input, integer ( kind = 4 ) MAXST, the maximum size available for STACK array;
!    should be about NPT to be safe, but MAX ( 10, 2*LOG2(NPT) ) usually
!    enough.
!
!    Input, real ( kind = 8 ) VCL(1:2,1:*), the coordinates of 2D vertices.
!
!    Input, integer ( kind = 4 ) IND(1:NPT), the indices in VCL of vertices to be
!    triangulated; vertices are inserted in order given by this array
!
!    Output, integer ( kind = 4 ) NTRI, the number of triangles in triangulation; equal to
!    2*NPT - NB - 2 where NB = number of boundary vertices.
!
!    Output, integer ( kind = 4 ) TIL(1:3,1:NTRI), the triangle incidence list; elements
!    are indices of VCL; vertices of triangles are in counterclockwise order.
!
!    Output, integer ( kind = 4 ) TNBR(1:3,1:NTRI), the triangle neighbor list;
!    positive elements are indices of TIL; negative elements are used for links
!    of counterclockwise linked list of boundary edges; LINK = -(3*I + J-1)
!    where I, J = triangle, edge index; TNBR(J,I) refers to
!    the neighbor along edge from vertex J to J+1 (mod 3).
!
!    Workspace, integer STACK(1:MAXST), the used for stack of triangles
!    for which circumcircle test must be made
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) maxst
  integer ( kind = 4 ) npt

  real    ( kind = 8 ) cmax
  integer ( kind = 4 ) e
  integer ( kind = 4 ) em1
  integer ( kind = 4 ) ep1
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i3
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
  integer ( kind = 4 ) m3
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
!
!  Determine initial triangle.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )

  m1 = ind(1)
  m2 = ind(2)

  do j = 1, 2
    cmax = max ( abs ( vcl(j,m1) ), abs ( vcl(j,m2) ) )
    if ( tol * cmax < abs ( vcl(j,m1) - vcl(j,m2) ) .and. tol < cmax ) then
      go to 20
    end if
  end do

  ierr = 224
  return
   20 continue

  i3 = 3
   30 continue

  if ( npt < i3 ) then
    ierr = 225
    return
  end if

  m = ind(i3)
  lr = lrline(vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1),vcl(1,m2), &
     vcl(2,m2),0.0D+00)

  if ( lr == 0 ) then
    i3 = i3 + 1
    go to 30
  end if

  if ( i3 /= 3 ) then
    ind(i3) = ind(3)
    ind(3) = m
  end if

  ntri = 1

  if ( lr == -1 ) then
    til(1,1) = m1
    til(2,1) = m2
  else
    til(1,1) = m2
    til(2,1) = m1
  end if

  til(3,1) = m
  tnbr(1,1) = -4
  tnbr(2,1) = -5
  tnbr(3,1) = -3

  if ( msglvl == 4 ) then
    write ( *,600) 1,vcl(1,m1),vcl(2,m1),vcl(1,m2),vcl(2,m2)
    write ( *,600) 1,vcl(1,m2),vcl(2,m2),vcl(1,m),vcl(2,m)
    write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
  end if
!
!  Insert vertices one at a time from anywhere, walk through
!  triangulation to determine location of new vertex, and apply
!  diagonal edge swaps until Delaunay triangulation of vertices
!  (so far) is obtained.
!
  top = 0

  do i = 4, npt

    if ( msglvl == 4 ) then
      write ( *,600) i
    end if

    m = ind(i)
    rtri = ntri
    call walkt2(vcl(1,m),vcl(2,m),ntri,vcl,til,tnbr,rtri,redg, ierr)

    if ( redg == 0 ) then

      m1 = til(1,rtri)
      m2 = til(2,rtri)
      m3 = til(3,rtri)
      til(3,rtri) = m

      if ( 0 < tnbr(1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m2
      til(2,ntri) = m3
      til(3,ntri) = m
      n = tnbr(2,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      ntri = ntri + 1
      til(1,ntri) = m3
      til(2,ntri) = m1
      til(3,ntri) = m
      n = tnbr(3,rtri)
      tnbr(1,ntri) = n

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      tnbr(2,rtri) = ntri - 1
      tnbr(3,rtri) = ntri
      tnbr(2,ntri-1) = ntri
      tnbr(3,ntri-1) = rtri
      tnbr(2,ntri) = rtri
      tnbr(3,ntri) = ntri - 1

      if ( tnbr(1,ntri-1) <= 0 ) then

        t = rtri
        e = 1

40      continue

        if ( 0 < tnbr(e,t) ) then

          t = tnbr(e,t)

          if ( til(1,t) == m2 ) then
            e = 3
          else if ( til(2,t) == m2 ) then
            e = 1
          else
            e = 2
          end if

          go to 40

        end if

        tnbr(e,t) = -3*ntri + 3

      end if

      if ( 0 <= tnbr(1,ntri) ) then

        t = ntri - 1
        e = 1

50      continue

        if ( 0 < tnbr(e,t) ) then

          t = tnbr(e,t)

          if ( til(1,t) == m3 ) then
            e = 3
          else if ( til(2,t) == m3 ) then
            e = 1
          else
            e = 2
          end if

          go to 50

        end if

        tnbr(e,t) = -3*ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m1),vcl(2,m1)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

    else if ( redg < 0 ) then

      redg = -redg
      ltri = 0

      call vbedg(vcl(1,m),vcl(2,m),vcl,til,tnbr,ltri,ledg,rtri, &
        redg)

      n = ntri + 1
      l = -tnbr(ledg,ltri)

60    continue

      t = l/3
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
        go to 100
      end if

      stack(top) = ntri

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m2),vcl(2,m2)
      end if

      if ( t /= rtri .or. e /= redg ) then
        go to 60
      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m), vcl(1,m1),vcl(2,m1)
      end if

      tnbr(ledg,ltri) = -3*n - 1
      tnbr(2,n) = -3*ntri - 2
      tnbr(3,ntri) = -l

    else if ( redg <= 3 ) then

      m1 = til(redg,rtri)

      if ( redg == 1 ) then
        e = 2
        ep1 = 3
      else if ( redg == 2 ) then
        e = 3
        ep1 = 1
      else
        e = 1
        ep1 = 2
      end if

      m2 = til(e,rtri)
      til(e,rtri) = m
      m3 = til(ep1,rtri)

      if ( 0 < tnbr(ep1,rtri) ) then
        top = 1
        stack(top) = rtri
      end if

      ntri = ntri + 1
      til(1,ntri) = m
      til(2,ntri) = m2
      til(3,ntri) = m3
      n = tnbr(e,rtri)
      tnbr(2,ntri) = n
      tnbr(3,ntri) = rtri
      tnbr(e,rtri) = ntri

      if ( 0 < n ) then

        if ( tnbr(1,n) == rtri ) then
          tnbr(1,n) = ntri
        else if ( tnbr(2,n) == rtri ) then
          tnbr(2,n) = ntri
        else
          tnbr(3,n) = ntri
        end if

        top = top + 1
        stack(top) = ntri

      end if

      if ( msglvl == 4 ) then
        write ( *,600) 1,vcl(1,m),vcl(2,m),vcl(1,m3),vcl(2,m3)
      end if

      ltri = tnbr(redg,rtri)

      if ( ltri <= 0 ) then

        tnbr(1,ntri) = ltri
        tnbr(redg,rtri) = -3*ntri
        if ( tnbr(2,ntri) <= 0) tnbr(1,ntri) = -3*ntri - 1

      else

        tnbr(1,ntri) = ntri + 1
        tnbr(redg,rtri) = ltri

        if ( til(1,ltri) == m2 ) then
          ledg = 1
          em1 = 2
          e = 3
        else if ( til(2,ltri) == m2 ) then
          ledg = 2
          em1 = 3
          e = 1
        else
          ledg = 3
          em1 = 1
          e = 2
        end if

        til(ledg,ltri) = m
        m3 = til(e,ltri)

        if ( 0 < tnbr(em1,ltri) ) then
          top = top + 1
          stack(top) = ltri
        end if

        ntri = ntri + 1
        til(1,ntri) = m2
        til(2,ntri) = m
        til(3,ntri) = m3
        tnbr(1,ntri) = ntri - 1
        tnbr(2,ntri) = ltri
        n = tnbr(e,ltri)
        tnbr(3,ntri) = n
        tnbr(e,ltri) = ntri

        if ( 0 < n ) then

          if ( tnbr(1,n) == ltri ) then
            tnbr(1,n) = ntri
          else if ( tnbr(2,n) == ltri ) then
            tnbr(2,n) = ntri
          else
            tnbr(3,n) = ntri
          end if

          top = top + 1
          stack(top) = ntri

        end if

        if ( msglvl == 4 ) then
          write ( *,600) 1,vcl(1,m),vcl(2,m), vcl(1,m3),vcl(2,m3)
        end if

        if ( tnbr(2,ntri-1) <= 0 ) then

          t = ntri
          e = 3

70        continue

          if ( 0 < tnbr(e,t) ) then

            t = tnbr(e,t)

            if ( til(1,t) == m2 ) then
              e = 3
            else if ( til(2,t) == m2 ) then
              e = 1
            else
              e = 2
            end if

            go to 70

          end if

          tnbr(e,t) = -3*ntri + 2

        end if

        if ( tnbr(3,ntri) <= 0 ) then

          t = ltri

          if ( ledg <= 2 ) then
            e = ledg + 1
          else
            e = 1
          end if

80        continue

          if ( 0 < tnbr(e,t) ) then

            t = tnbr(e,t)

            if ( til(1,t) == m3 ) then
              e = 3
            else if ( til(2,t) == m3 ) then
              e = 1
            else
              e = 2
            end if

            go to 80

          end if

          tnbr(e,t) = -3*ntri - 2

        end if

      end if

    else
      ierr = 224
      go to 100
    end if

    call swapec ( m, top, maxst, 0, 0, vcl, til, tnbr, stack, ierr )

    if ( ierr /= 0 ) then
      exit
    end if

  end do

100 continue

  if ( i3 /= 3 ) then
    t = ind(i3)
    ind(i3) = ind(3)
    ind(3) = t
  end if

  if ( msglvl == 4 ) then
    write ( *,600) npt+1
  end if

  600 format (1x,i7,4f15.7)

  return
end
