subroutine vprght ( xc, yc, ivis, ierr )

!*****************************************************************************80
!
!! VPRGHT is called by VISPOL for the RIGHT operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the RIGHT
!    operation (OPER = 2).
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
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:*)
  real    ( kind = 8 ) xe
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xw
  real    ( kind = 8 ) yc(0:*)
  real    ( kind = 8 ) ye
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, EYE-S(TOP)-V(CUR) is a right
!  turn, EYE-S(TOP-1)-S(TOP) is a left turn, TOP < CUR, S(TOP) =
!  V(CUR-1) and S(TOP-1)-S(TOP)-V(CUR) is a left turn or S(TOP) is
!  not on edge V(CUR-1)-V(CUR) and V(CUR-1)-V(CUR) intersects
!  EYE-S(TOP).
!  Pop points from stack. If BEYE, it is not possible for
!  (XC(CUR),YC(CUR)) to be identical to any stack points.
!
  ierr = 0

10 continue

  case = 0
  j = top

20 continue

  if ( abs(ivis(j)) <= nv ) then

    lr = lrline(xc(cur),yc(cur),xe,ye,xc(j-1),yc(j-1),0.0D+00)

    if ( lr == -1 ) then

      case = 1

    else if ( lr == 0 ) then

      if ( abs(ivis(j-1)) <= nv ) then

        j = j - 1
        case = 2

      else if ( (xc(j-1) - xe)**2 + (yc(j-1) - ye)**2 <= &
                (xc(j-2) - xe)**2 + (yc(j-2) - ye)**2 ) then

        j = j - 2
        case = 2

      else

        case = -1

      end if

    end if

  else if ( case == -1 ) then

    if ( (xc(cur) - xe)**2 + (yc(cur) - ye)**2 <= &
         (xc(j-1) - xe)**2 + (yc(j-1) - ye)**2 ) then
      j = j - 1
      case = 2
    else
      xw = xc(cur)
      yw = yc(cur)
      case = 3
    end if

  else

    call xedge(0,xc(cur-1),yc(cur-1),xc(cur),yc(cur), &
           xc(j-1),yc(j-1),xc(j),yc(j),xw,yw,intsct)

    if ( intsct ) then
      case = 3
    end if

  end if

  if ( 0 < case ) go to 30
  j = j - 1
  if ( 1 <= j ) go to 20
!
!  Error from no more edges in stack.
!
  ierr = 206
  return
!
!  Process next edge.
!
30 continue

  if ( case == 3 ) then

    oper = 6
    top = j - 1

  else

    top = j
    xw = xc(cur-1)
    yw = yc(cur-1)

    if ( case == 1 ) then
      call xedge(1,xe,ye,xc(cur),yc(cur),xc(top-1),yc(top-1), &
           xc(top),yc(top),xu,yu,intsct)
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -abs(ivis(top))
    end if

    lr = lrline(xc(cur+1),yc(cur+1),xe,ye,xc(cur),yc(cur),0.0D+00)

    if ( lr == 1 ) then

      cur = cur + 1

    else

      j = cur + 1
      lr1 = lrline(xc(j),yc(j),xw,yw,xc(cur),yc(cur),0.0D+00)

      if ( lr1 == -1 ) then

        oper = 5
        cur = j

      else

        if ( lr == -1 ) then
          lr2 = -1
          go to 50
        end if

40      continue

        j = j + 1
        lr2 = lrline(xc(j),yc(j),xe,ye,xc(cur),yc(cur),0.0D+00)
        if ( lr2 == 0) go to 40

50      continue

        if ( lr2 == -1 ) then
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

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 2) go to 10

  return
end
