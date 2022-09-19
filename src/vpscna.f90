subroutine vpscna ( xc, yc, ivis, ierr )

!*****************************************************************************80
!
!! VPSCNA is called by VISPOL for the SCANA operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the SCANA
!    operation (OPER = 3).
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
!    Input/output, real ( kind = 8 ) XC, YC, see comments in VISPOL.
!
!    Input/output, integer ( kind = 4 ) IVIS, see comments in VISPOL.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  logical              beye
  integer ( kind = 4 ) case
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierr
  logical              intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lr3
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:*)
  real    ( kind = 8 ) xe
  real    ( kind = 8 ) xw
  real    ( kind = 8 ) yc(0:*)
  real    ( kind = 8 ) ye
  real    ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn or forward move, S(TOP) =
!  V(CUR-1) or EYE-S(TOP)-V(CUR-1) is a forward move and TOP = 0,
!  TOP < CUR; S(TOP-1)-S(TOP)-V(CUR) is a right turn if 1 <= TOP
!  or EYE-S(TOP)-V(CUR) is a right turn if TOP = 0.
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  vertex but any edge incident on this vertex encountered during
!  scan must be invisible from (XE,YE).
!
  ierr = 0
  k = cur

10 continue

  if ( xc(k+1) == xc(top) .and. yc(k+1) == yc(top) ) then

    k = k + 2

  else

    call xedge(1,xe,ye,xc(top),yc(top),xc(k),yc(k),xc(k+1), &
      yc(k+1),xw,yw,intsct)

    if ( intsct ) then

      lr = lrline(xc(k+1),yc(k+1),xe,ye,xc(k),yc(k),0.0D+00)

      if ( lr == 1 ) then

        if ( (xw - xe)**2 + (yw - ye)**2 <= &
             (xc(top) - xe)**2 + (yc(top) - ye)**2 ) then

          if ( 0 < top ) then
            case = 1
            go to 20
          end if

        else

          lr1 = lrline(xc(k),yc(k),xe,ye,xc(top),yc(top),0.0D+00)

          if ( lr1 == -1 ) then
            case = 2
            go to 20
          end if

        end if

      else

        lr1 = lrline(xc(k+1),yc(k+1),xe,ye,xc(top),yc(top), 0.0D+00)

        if ( lr1 == -1 ) then
          case = 3
          go to 20
        end if

      end if

    end if

    k = k + 1

  end if

  if ( k < nv) go to 10
!
!  Error from unsuccessful scan.
!
  ierr = 207
  return
!
!  Process current edge.
!
20 continue

  if ( case == 3 ) then

    oper = 1
    cur = k + 1
    lr = lrline(xc(k),yc(k),xe,ye,xc(top),yc(top),0.0D+00)
    top = top + 1

    if ( lr == 0 ) then
      xc(top) = xc(k)
      yc(top) = yc(k)
      ivis(top) = k + nv
    else
      xc(top) = xw
      yc(top) = yw
      ivis(top) = -(k + 1 + nv)
    end if

    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else if ( case == 1 ) then

    cur = k + 1
    lr = lrline(xc(cur),yc(cur),xe,ye,xc(top),yc(top),0.0D+00)

    if ( lr == 1 ) then

      oper = 2

    else

      j = cur + 1
      lr1 = lrline(xc(j),yc(j),xe,ye,xc(cur),yc(cur),0.0D+00)
      lr2 = lrline(xc(j),yc(j),xc(k),yc(k),xc(cur),yc(cur),0.0D+00)

      if ( lr1 <= 0 .and. lr2 == -1 ) then

        oper = 5
        xw = xc(k)
        yw = yc(k)
        cur = j

      else

        if ( lr1 /= 0 ) then
          lr3 = lr1
          go to 40
        end if

30      continue

        j = j + 1
        lr3 = lrline(xc(j),yc(j),xe,ye,xc(cur),yc(cur),0.0D+00)
        if ( lr3 == 0) go to 30

40      continue

        if ( lr3 == 1 ) then
          oper = 2
        else
          oper = 1
          top = top + 1
          xc(top) = xc(j-1)
          yc(top) = yc(j-1)
          ivis(top) = j - 1 + nv
          top = top + 1
          xc(top) = xc(j)
          yc(top) = yc(j)
          ivis(top) = j
        end if

        cur = j

      end if

    end if

  else

    oper = 6
    cur = k + 1
    lr = lrline(xc(cur),yc(cur),xe,ye,xc(top),yc(top),0.0D+00)

    if ( lr == 0 ) then
      xw = xc(cur)
      yw = yc(cur)
    end if

  end if

  return
end
