subroutine vispol ( xeye, yeye, nvrt, xc, yc, nvis, ivis, ierr )

!*****************************************************************************80
!
!! VISPOL computes the visibility polygon from an eyepoint.
!
!  Discussion:
!
!    This routine computes the visibility polygon VP from an eyepoint in
!    the interior or blocked exterior of a simple polygon P or
!    on the boundary of a simply connected polygonal region P.
!    In the latter case, the interior angles at all vertices must
!    be strictly between 0 and 2*PI.
!
!    On input, XC and YC contain vertex coordinates of P. During
!    the algorithm, part of XC, YC is used as a stack, which, on
!    output, contains the vertex coordinates of VP.  The stack
!    vertices overwrite the input vertices as the input vertices
!    are scanned.  Elements of IVIS are set when vertices are added
!    to the stack; these values may have +NV or -NV added to them
!    to indicate that stack point has same angle as previous one.
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
!    Barry Joe and R. B. Simpson,
!    BIT
!    Volume 27, 1987, pages 458-473.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XEYE, YEYE, the coordinates of eyepoint; must
!    be a simple vertex if it lies on the boundary (i.e. occurs only once).
!
!    Input, integer ( kind = 4 ) NVRT, the upper subscript of XC, YC (approximate number
!    of vertices).
!
!    Input/output, real ( kind = 8 ) XC(0:NVRT), YC(0:NVRT).  On input, if
!    eyepoint is interior or blocked exterior, then arrays contain
!    coordinates in counterclockwise or CW order, respectively, with
!    (XC(0),YC(0)) = (XC(NVRT), YC(NVRT)); (XC(0),YC(0)) is a vertex visible
!    from (XEYE,YEYE), e.g. as computed by routine ROTIPG.
!    If eyepoint is a vertex of P then arrays contain
!    coordinates in counterclockwise order; (XC(0),YC(0)) is successor
!    vertex of (XEYE,YEYE); (XC(NVRT),YC(NVRT)) is
!    predecessor vertex of (XEYE,YEYE).
!    On output, vertices of VP in counterclockwise order;
!    if eyepoint is interior or blocked exterior then
!    (XC(0),YC(0)) = (XC(NVIS),YC(NVIS)), else (XC(0),YC(0))
!    and (XC(NVIS),YC(NVIS)) are the successor and
!    predecessor vertices of (XEYE,YEYE) in VP.
!
!    Output, integer ( kind = 4 ) NVIS, the upper subscript of XC, YC on output (approximate
!    number of vertices of VP); NVIS <= NVRT.
!
!    Output, integer ( kind = 4 ) IVIS(0:NVIS), contains information about the vertices
!    of VP with respect to the vertices of P; IVIS(I) = K if ( XC(I),YC(I))
!    is the vertex of index K in the input polygon; IVIS(I)
!    = -K if ( XC(I),YC(I)) is on the interior of the edge
!    joining vertices of index K-1 and K in input polygon.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) nvrt

  logical              beye
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ivis(0:nvrt)
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) nvis
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:nvrt)
  real    ( kind = 8 ) xe
  real    ( kind = 8 ) xeye
  real    ( kind = 8 ) xw
  real    ( kind = 8 ) yc(0:nvrt)
  real    ( kind = 8 ) ye
  real    ( kind = 8 ) yeye
  real    ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  Variables in common block GVPVAR:
!
!        NV - NVRT
!        OPER - operation code 1 to 7 for LEFT, RIGHT, SCANA, SCANB,
!              SCANC, SCAND, FINISH
!        CUR - index of current vertex of P in XC, YC arrays
!        TOP - index of top vertex of stack in XC, YC arrays
!              (TOP <= CUR is always satisfied)
!        XE,YE - XEYE,YEYE
!        XW,YW - coordinates of point on last or second-last edge
!              processed (needed for routines VPSCNB, VPSCNC, VPSCND)
!        BEYE - TRUE iff eyepoint is on boundary
!
  ierr = 0
  beye = xc(0) /= xc(nvrt) .or. yc(0) /= yc(nvrt)
  nv = nvrt
  xe = xeye
  ye = yeye
  ivis(0) = 0
  cur = 1

  if ( beye ) then

    do

      lr = lrline ( xc(nv-1), yc(nv-1), xe, ye, xc(nv), yc(nv), 0.0D+00 )

      if ( lr /= 0 ) then
        exit
      end if

      nv = nv - 1

    end do

  end if

20 continue

  lr = lrline(xc(cur),yc(cur),xe,ye,xc(0),yc(0),0.0D+00)

  if ( lr == 0 ) then
    cur = cur + 1
    go to 20
  end if

  if ( lr == -1 ) then

    oper = 1

    if ( cur == 1 ) then

      top = 1
      ivis(1) = cur

    else if ( beye ) then

      top = 1
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
      xc(1) = xc(cur)
      yc(1) = yc(cur)
      ivis(1) = cur

    else

      top = 2
      xc(1) = xc(cur-1)
      yc(1) = yc(cur-1)
      ivis(1) = cur - 1 + nv
      xc(2) = xc(cur)
      yc(2) = yc(cur)
      ivis(2) = cur

    end if

  else

    oper = 3
    top = 0

    if ( beye .and. 1 < cur ) then
      xc(0) = xc(cur-1)
      yc(0) = yc(cur-1)
      ivis(0) = cur - 1
    end if

  end if
!
!  Angular displacement of stack points are in nondecreasing order,
!  with at most two consecutive points having the same displacement.
!
30 continue

  if ( oper == 1 ) then
    call vpleft(xc,yc,ivis)
  else if ( oper == 2 ) then
    call vprght(xc,yc,ivis, ierr )
  else if ( oper == 3 ) then
    call vpscna(xc,yc,ivis, ierr )
  else if ( oper == 4 ) then
    call vpscnb(xc,yc,ivis, ierr )
  else if ( oper == 5 ) then
    call vpscnc(xc,yc,ivis, ierr )
  else
    call vpscnd(xc,yc,ivis, ierr )
  end if

  if ( ierr /= 0 ) then
    nvis = top
    return
  end if

  if ( oper <= 6 ) then
    go to 30
  end if
!
!  Add or subtract NV from those IVIS values which are used to
!  indicate that stack point has same angle as previous one.
!
  do i = 1, top
    if ( nv < ivis(i) ) then
      ivis(i) = ivis(i) - nv
    else if ( ivis(i) < -nv ) then
      ivis(i) = ivis(i) + nv
    end if
  end do

  nvis = top

  return
end
