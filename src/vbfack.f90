subroutine vbfack ( k, pt, ctr, vcl, vm, bf, fc, topv, topnv, ind, mat, vec )

!*****************************************************************************80
!
!! VBFACK determines the boundary faces of a KD triangulation.
!
!  Discussion:
!
!    This routine determines boundary faces of a K-D triangulation visible
!    from point PT, given a starting visible boundary face.
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
!    Input, integer ( kind = 4 ) K, the dimension of triangulation.
!
!    Input, real ( kind = 8 ) PT(1:K), the K-D point.
!
!    Input, real ( kind = 8 ) CTR(1:K), the K-D point in interior of
!    triangulation.
!
!    Input, real ( kind = 8 ) VCL(1:K,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:*), the indices of vertices of VCL being triangulated.
!
!    Input, integer ( kind = 4 ) BF(1:K,1:*), the  array of boundary face records;
!    see DTRISK.
!
!    Input/output, integer ( kind = 4 ) FC(1:K+4,1:*), the array of face records; see routine
!    DTRISK; row K+4 is used for links of 3 stacks in this routine.
!    On output, FC(K+4,*) - gets updated; only stack of visible boundary faces
!    remains at end of routine.
!
!    Input/output, integer ( kind = 4 ) TOPV.  On input, index of FC of visible boundary
!    face.  On output, index of top of stack of visible boundary faces.
!
!    Input, integer ( kind = 4 ) TOPNV, the index of top of stack of boundary faces
!    already found to be not visible from PT, or 0 for empty stack.
!
!    Workspace, integer IND(1:K), the indices in VCL of K-D vertices.
!
!    Workspace, real ( kind = 8 ) MAT(1:K-1,1:K), the matrix used for
!    solving system of homogeneous linear equations.
!
!    Workspace, real ( kind = 8 ) VEC(1:K), the vector used for
!    hyperplane normal.
!
  implicit none

  integer ( kind = 4 ) k

  integer ( kind = 4 ) bf(k,*)
  real    ( kind = 8 ) ctr(k)
  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) ind(k)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp2
  integer ( kind = 4 ) kp4
  integer ( kind = 4 ) l
  real    ( kind = 8 ) mat(k-1,k)
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) opsidk
  real    ( kind = 8 ) pt(k)
  integer ( kind = 4 ) ptr
  integer ( kind = 4 ) topn
  integer ( kind = 4 ) topnv
  integer ( kind = 4 ) topt
  integer ( kind = 4 ) topv
  real    ( kind = 8 ) vcl(k,*)
  real    ( kind = 8 ) vec(k)
  integer ( kind = 4 ) vm(*)
!
!  TOPN is index of top of stack of non-visible boundary faces.
!  TOPT is index of top of stack of boundary faces to be tested.
!
  kp2 = k + 2
  kp4 = k + 4
  topn = topnv
  topt = 0
  fc(kp4,topv) = 0
  l = -fc(kp2,topv)

  do j = 1,k
    nbr = bf(j,l)
    if ( fc(kp4,nbr) == -1 ) then
      fc(kp4,nbr) = topt
      topt = nbr
    end if
  end do

  do while ( topt /= 0 )

    ptr = topt
    topt = fc(kp4,ptr)
    ind(1:k) = vm(fc(1:k,ptr))

    if ( opsidk(k,ind,vcl,.false.,ctr,pt,mat,vec) == 1 ) then

      fc(kp4,ptr) = topv
      topv = ptr
      l = -fc(kp2,ptr)

      do j = 1,k
        nbr = bf(j,l)
        if ( fc(kp4,nbr) == -1 ) then
          fc(kp4,nbr) = topt
          topt = nbr
        end if
      end do

    else

      fc(kp4,ptr) = topn
      topn = ptr

    end if

  end do
!
!  For boundary faces not visible from PT, set FC(KP4,*) = -1.
!
  do while (topn /= 0)
    ptr = topn
    topn = fc(kp4,ptr)
    fc(kp4,ptr) = -1
  end do

  return
end
