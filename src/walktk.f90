subroutine walktk ( k, pt, n, p, nsmplx, vcl, vm, fc, ht, ifac, ivrt, indl, &
  indv, ipvt, alpha, mat, ierr )

!*****************************************************************************80
!
!! WALKTH finds the Delaunay simplex containing a point by "walking".
!
!  Discussion:
!
!    This routine walks through neighboring simplices of a K-D (Delaunay)
!    triangulation until a simplex is found containing point PT
!    or PT is found to be outside the convex hull. Search is
!    guaranteed to terminate for a unique Delaunay triangulation,
!    else a cycle may occur.
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
!    Input, integer ( kind = 4 ) N, the upper bound on vertex indices or size of VM.
!
!    Input, integer ( kind = 4 ) P, the size of hash table.
!
!    Input, integer ( kind = 4 ) NSMPLX, the number of simplices in triangulation; used
!    to detect cycle.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) VM(1:N), the vertex mapping list (maps from local
!    indices used in FC to indices of VCL).
!
!    Input, integer ( kind = 4 ) FC(1:K+4,1:*), array of face records; see routine DTRISK.
!
!    Input, integer ( kind = 4 ) HT(0:P-1), the hash table using direct chaining.
!
!    Input/output, integer ( kind = 4 ) IFAC.  On input, index of face of FC to begin
!    search at.  On ouput, index of FC indicating simplex or face containing PT.
!
!    Input/output, integer ( kind = 4 ) IVRT.  On input, K+1 or K+2 to indicate simplex to
!    begin search at.  On output, K+1 or K+2 to indicate that FC(IVRT,IFAC) is
!    (K+1)st vertex of simplex containing PT in its interior; 0 if PT
!    is outside convex hull on other side of face FC(*,IFAC);
!    K if PT lies in interior of face FC(*,IFAC); 1 to K-1 if
!    PT lies in interior of facet of FC(*,IFAC) of dim IVRT-1.
!
!    Output, integer ( kind = 4 ) INDV(1:K-1), if 1 <= IVRT <= K-1 then first IVRT elements
!    are local vertex indices in increasing order of (IVRT-1)-
!    facet containing PT in its interior.
!
!    Workspace, integer INDL(1:K), the local vertex indices of K-D vertices.
!
!    Workspace, integer INDV(1:K+1), the indices in VCL of K-D vertices.
!
!    Workspace, integer IPVT(1:K-1), the pivot indices.
!
!    Workspace, real ( kind = 8 ) ALPHA(1:K+1), the barycentric coordinates.
!
!    Workspace, real ( kind = 8 ) MAT(1:K,1:K), the matrix used for solving
!    system of linear equations.
!
!    Output, integer ( kind = 4 ) IERR, error flag, which is zero unless an error occurred.
!
  implicit none

  integer ( kind = 4 ) k
  integer ( kind = 4 ) n
  integer ( kind = 4 ) p

  integer ( kind = 4 ) a
  real    ( kind = 8 ) alpha(k+1)
  integer ( kind = 4 ) cnt
  integer ( kind = 4 ) d
  logical              degen
  integer ( kind = 4 ) fc(k+4,*)
  integer ( kind = 4 ) ht(0:p-1)
  integer ( kind = 4 ) htsrck
  integer ( kind = 4 ) i
  integer ( kind = 4 ) ierr
  integer ( kind = 4 ) ifac
  integer ( kind = 4 ) indl(k)
  integer ( kind = 4 ) indv(k+1)
  integer ( kind = 4 ) ineg
  integer ( kind = 4 ) ipvt(k-1)
  integer ( kind = 4 ) ivrt
  integer ( kind = 4 ) izero
  integer ( kind = 4 ) j
  integer ( kind = 4 ) kp1
  integer ( kind = 4 ) kp2
  real    ( kind = 8 ) mat(k,k)
  integer ( kind = 4 ) npos
  integer ( kind = 4 ) nsmplx
  integer ( kind = 4 ) nzero
  real    ( kind = 8 ) pt(k)
  real    ( kind = 8 ) tol
  real    ( kind = 8 ) vcl(k,*)
  integer ( kind = 4 ) vm(n)

  ierr = 0
  tol = 100.0D+00 * epsilon ( tol )
  kp1 = k + 1
  kp2 = k + 2
  cnt = 0

10 continue

  if ( fc(ivrt,ifac) <= 0 ) then
    ivrt = 0
    return
  end if

  cnt = cnt + 1

  if ( nsmplx < cnt ) then
    ierr = 407
    return
  end if

  indl(1:k) = fc(1:k,ifac)
  indv(1:k) = vm(indl(1:k))

  d = fc(ivrt,ifac)
  indv(kp1) = vm(d)
  call baryck(k,indv,vcl,pt,alpha,degen,mat,ipvt)

  if ( degen ) then
    ierr = 401
    return
  end if

  npos = 0
  nzero = 0
  izero = 0
  ineg = 0

  do i = 1, kp1

    if ( tol < alpha(i) ) then
      npos = npos + 1
    else if ( alpha(i) < -tol ) then
      ineg = i
    else
      nzero = nzero + 1
      izero = i
    end if

  end do

  if ( npos == kp1 ) then

    return

  else if ( ineg == kp1 ) then

    ivrt = kp1 + kp2 - ivrt

  else if ( 0 < ineg ) then

    a = indl(ineg)
    indl(ineg) = d
    ifac = htsrck(k,indl,n,p,fc,ht)

    if ( ifac <= 0 ) then
      ierr = 400
      return
    end if

    if ( fc(kp1,ifac) == a ) then
      ivrt = kp2
    else
      ivrt = kp1
    end if

  else if ( nzero == 1 ) then

    ivrt = k

    if ( izero < kp1 ) then
      indl(izero) = d
      ifac = htsrck(k,indl,n,p,fc,ht)
      if ( ifac <= 0) ierr = 400
    end if

    return

  else

    ivrt = npos
    j = 0

    do i = 1, k
      if ( tol < alpha(i) ) then
        j = j + 1
        indv(j) = indl(i)
      end if
    end do

    if ( izero < kp1 ) then

      j = npos

50    continue

      if ( 1 < j ) then

        if ( d < indv(j-1) ) then
          indv(j) = indv(j-1)
          j = j - 1
          go to 50
        end if

      end if

      indv(j) = d
      indl(izero) = d
      ifac = htsrck(k,indl,n,p,fc,ht)

      if ( ifac <= 0 ) then
        ierr = 400
        return
      end if

    end if

    return

  end if

  go to 10

end
