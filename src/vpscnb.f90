subroutine vpscnb ( xc, yc, ivis, ierr )

!*****************************************************************************80
!
!! VPSCNB is called by VISPOL for the SCANB operation.
!
!  Discussion:
!
!    This routine is called by routine VISPOL for the SCANB
!    operation (OPER = 4).
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
  integer ( kind = 4 ) lrline
  integer ( kind = 4 ) nv
  integer ( kind = 4 ) oper
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) tolabs
  integer ( kind = 4 ) top
  real    ( kind = 8 ) xc(0:*)
  real    ( kind = 8 ) xe
  real    ( kind = 8 ) xu
  real    ( kind = 8 ) xw
  real    ( kind = 8 ) yc(0:*)
  real    ( kind = 8 ) ye
  real    ( kind = 8 ) yw
  real    ( kind = 8 ) yu

  common /gvpvar/ nv,oper,cur,top,xe,ye,xw,yw,beye
  save /gvpvar/
!
!  EYE-V(CUR-1)-V(CUR) is a left turn, S(TOP) = V(CUR) or S(TOP) is
!  on interior of edge V(CUR-1)-V(CUR), TOP <= CUR, S(TOP) has
!  angular displacement of 2*PI or interior angle at boundary eye.
!  (XW,YW) is the input version of (XC(CUR),YC(CUR)).
!  If BEYE, it is possible that (XC(TOP),YC(TOP)) is a non-simple
!  point but any edge containing this point encountered during scan
!  must be invisible from (XE,YE), except for 1 case where K = CUR.
!
  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  tolabs = tol*((xc(nv) - xc(top))**2 + (yc(nv) - yc(top))**2)
  k = cur

  if ( ivis(top) < 0 .or. k + 1 == nv ) then
    go to 10
  end if

  lr = lrline(xc(k+1),yc(k+1),xe,ye,xc(top),yc(top),0.0D+00)
  lr1 = lrline(xc(k+1),yc(k+1),xc(top-1),yc(top-1),xc(top),yc(top), 0.0D+00)

  if ( lr == 1 .and. lr1 == -1 ) then
    oper = 2
    cur = k + 1
    return
  else
    k = k + 1
  end if

10 continue

  if ( k + 1 == nv ) then

    oper = 7
    cur = nv
    top = top + 1
    xc(top) = xc(nv)
    yc(top) = yc(nv)
    ivis(top) = nv
    return

  else

    if ( k == cur ) then
      call xedge(0,xc(nv),yc(nv),xc(top),yc(top),xw,yw, &
        xc(k+1),yc(k+1),xu,yu,intsct)
    else
      call xedge(0,xc(nv),yc(nv),xc(top),yc(top),xc(k),yc(k), &
        xc(k+1),yc(k+1),xu,yu,intsct)
    end if

    if ( intsct ) then

      if ( (xc(top) - xu)**2 + (yc(top) - yu)**2 <= tolabs ) then
        go to 20
      end if

      lr = lrline(xc(k+1),yc(k+1),xe,ye,xc(nv),yc(nv),0.0D+00)

      if ( lr == 1 ) then
        oper = 2
        cur = k + 1
        return
      end if

    end if

20  continue

    k = k + 1

  end if

  if ( k < nv ) then
    go to 10
  end if
!
!  Error from unsuccessful scan.
!
  ierr = 208

  return
end
