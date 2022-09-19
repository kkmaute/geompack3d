subroutine vpscnc ( xc, yc, ivis, ierr )

!*****************************************************************************80
!
!! VPSCNC is called by VISPOL for the SCANC operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the SCANC
!    operation (OPER = 5).
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
  integer ( kind = 4 ) j
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
!  EYE-V(CUR-1)-V(CUR) is a left turn or forward move, EYE-V(CUR-2)-
!  V(CUR-1) is a right turn, V(CUR-2)-V(CUR-1)-V(CUR) is a left turn,
!  TOP < CUR-1, W = V(CUR-2), S(TOP) is not on V(CUR-1)-V(CUR), EYE-
!  S(TOP)-V(CUR-1) is a backward move, EYE-S(TOP-1)-S(TOP) is a left
!  turn. If BEYE, it is possible that V(CUR-1) is a non-simple point,
!  but intersection at (XC(TOP),YC(TOP)) cannot occur.
!
  ierr = 0
  xp = xc(cur-1)
  yp = yc(cur-1)
  k = cur

10 continue

  if ( xc(k+1) == xp .and. yc(k+1) == yp ) then

    go to 40

  else if ( xc(k) == xp .and. yc(k) == yp ) then

    j = k + 1
    lr = lrline(xc(j),yc(j),xe,ye,xp,yp,0.0D+00)
    lr1 = lrline(xc(j),yc(j),xw,yw,xp,yp,0.0D+00)

    if ( lr <= 0 .and. lr1 == -1 ) go to 40

    if ( lr /= 0 ) then
      lr2 = lr
      go to 30
    end if

20  continue

    j = j + 1
    lr2 = lrline(xc(j),yc(j),xe,ye,xp,yp,0.0D+00)

    if ( lr2 == 0) go to 20

30  continue

    if ( lr2 == 1 ) then
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
    return

  else

    call xedge(0,xp,yp,xc(top),yc(top),xc(k),yc(k),xc(k+1), &
      yc(k+1),xu,yu,intsct)

    if ( intsct ) then

      lr = lrline(xc(k+1),yc(k+1),xe,ye,xp,yp,0.0D+00)

      if ( lr == 1 ) then
        oper = 2
        cur = k + 1
        return
      end if

    end if

  end if

40 continue

  k = k + 1
  if ( k < nv) go to 10
!
!  Error from unsuccessful scan.
!
  ierr = 209

  return
end
