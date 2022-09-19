subroutine vpscnd ( xc, yc, ivis, ierr )

!*****************************************************************************80
!
!! VPSCND is called by VISPOL for the SCAND operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the SCAND
!    operation (OPER = 6).
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
  integer ( kind = 4 ) cur
  integer ( kind = 4 ) ierr
  logical              intsct
  integer ( kind = 4 ) ivis(0:*)
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lr
  integer ( kind = 4 ) lr1
  integer ( kind = 4 ) lr2
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:*)
  real    ( kind = 8 ) xe
  real    ( kind = 8 ) xp
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xw
  real    ( kind = 8 ) yc(0:*)
  real    ( kind = 8 ) ye
  real    ( kind = 8 ) yp
  real    ( kind = 8 ) yu
  real    ( kind = 8 ) yw

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a right turn, S(TOP) is a V vertex not on
!  V(CUR-1)-V(CUR), TOP < CUR, W is intersection of V(CUR-1)-V(CUR)
!  and ray EYE-S(TOP), EYE-S(TOP)-W is a forward move, and
!  EYE-S(TOP-1)-S(TOP) is a left turn if 1 <= TOP.
!  If BEYE, it is possible that (XW,YW) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierr = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  call xedge(0,xw,yw,xc(top),yc(top),xc(k),yc(k),xc(k+1), &
    yc(k+1),xu,yu,intsct)

  if ( intsct ) then

    lr = lrline(xc(k+1),yc(k+1),xe,ye,xc(k),yc(k),0.0D+00)
    lr1 = lrline(xc(k+1),yc(k+1),xe,ye,xc(top),yc(top),0.0D+00)

    if ( lr == -1 .and. lr1 == -1 ) then

      if ( xc(k) /= xw .or. yc(k) /= yw ) go to 20

      lr2 = lrline(xc(k+1),yc(k+1),xp,yp,xw,yw,0.0D+00)
      if ( lr2 == -1) go to 30

20    continue

      oper = 1
      cur = k + 1
      lr2 = lrline(xc(k),yc(k),xe,ye,xc(top),yc(top),0.0D+00)
      top = top + 1

      if ( lr2 == 0 ) then
        xc(top) = xc(k)
        yc(top) = yc(k)
        ivis(top) = k + nv
      else
        xc(top) = xu
        yc(top) = yu
        ivis(top) = -(k + 1 + nv)
      end if

      top = top + 1
      xc(top) = xc(cur)
      yc(top) = yc(cur)
      ivis(top) = cur
      return

    end if

  end if

30 continue

  k = k + 1
  if ( k < nv) go to 10
!
!  Error from unsuccessful scan.
!
  ierr = 210

  return
end
