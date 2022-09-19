subroutine vpleft ( xc, yc, ivis )

!*****************************************************************************80
!
!! VPLEFT is called by VISPOL for the LEFT operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the LEFT
!    operation (OPER = 1).
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
!    Input/output, real ( kind = 8 ) XC, YC, see comments in routine VISPOL.
!
!    Input/output, integer ( kind = 4 ) IVIS(0:*), see comments in routine VISPOL.
!
  implicit none

  logical              beye
  integer ( kind = 4 ) cur
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
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR), TOP <= CUR,
!  S(TOP-1) = V(CUR-1) or on interior of edge V(CUR-1)-V(CUR).
!
10 continue

  if ( cur == nv ) then
    oper = 7
    return
  end if

  if ( ( .not. beye ) .and. top <= 2 ) go to 20
!
!  Check if angular displacement of stack chain greater than or equal
!  to 2*PI or interior angle at boundary viewpoint.
!
  call xedge(1,xe,ye,xc(nv),yc(nv),xc(top-1),yc(top-1),xc(top), &
     yc(top),xu,yu,intsct)

  if ( intsct ) then

    oper = 4
    xw = xc(cur)
    yw = yc(cur)
    lr = lrline(xc(top),yc(top),xe,ye,xc(nv),yc(nv),0.0D+00)

    if ( lr == -1 ) then
      xc(top) = xu
      yc(top) = yu
      ivis(top) = -cur
    end if

    return

  end if
!
!  Process next edge.
!
20 continue

  lr = lrline(xc(cur+1),yc(cur+1),xe,ye,xc(cur),yc(cur),0.0D+00)

  if ( lr == -1 ) then

    cur = cur + 1
    top = top + 1
    xc(top) = xc(cur)
    yc(top) = yc(cur)
    ivis(top) = cur

  else

    j = cur + 1
    lr1 = lrline(xc(j),yc(j),xc(top-1),yc(top-1),xc(cur),yc(cur), 0.0D+00)

    if ( lr1 == 1 ) then

      oper = 3
      cur = j

    else

      if ( lr == 1 ) then
        lr2 = 1
        go to 40
      end if

30    continue

      j = j + 1
      lr2 = lrline(xc(j),yc(j),xe,ye,xc(cur),yc(cur),0.0D+00)

      if ( lr2 == 0 ) then
        go to 30
      end if

40    continue

      if ( lr2 == -1 ) then
        top = top + 1
        xc(top) = xc(j-1)
        yc(top) = yc(j-1)
        ivis(top) = j - 1 + nv
        top = top + 1
        xc(top) = xc(j)
        yc(top) = yc(j)
        ivis(top) = j
      else
        oper = 2
      end if

      cur = j

    end if

  end if
!
!  This test avoids extra subroutine calls.
!
  if ( oper == 1 ) then
    go to 10
  end if

  return
end
