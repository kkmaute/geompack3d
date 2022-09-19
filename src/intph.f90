subroutine intph ( hflag, umdf, headp, widp, nfcev, nedev, nvrev, listev, &
  infoev, ivrt, facval, edgval, vrtval, vcl, facep, fvl, pfl, cntr, mdfint, &
  mean, stdv, volp, nf, indf, meanf, stdvf )

!*****************************************************************************80
!
!! INTPH integrates a mesh distribution function over a polyhedron.
!
!  Discussion:
!
!    This routine computes the integral of MDF3(X,Y,Z) [heuristic mdf] or
!    UMDF(X,Y,Z) [user-supplied mdf] in convex polyhedron.
!
!    Parameters WIDP to VRTVAL are used only if HFLAG = TRUE.
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
!    Input, logical HFLAG, TRUE if heuristic mdf, FALSE if user-supplied mdf.
!
!    Input, external real ( kind = 8 ) UMDF(X,Y,Z), the user-supplied mdf
!    with d.p arguments.
!
!    Input, integer ( kind = 4 ) HEADP, the head pointer to face of PFL for convex
!    polyhedron.
!
!    Input, real ( kind = 8 ) WIDP, the width of original polyhedron of
!    decomposition.
!
!    Input, integer ( kind = 4 ) NFCEV, NEDEV, NVREV, LISTEV(1:NFCEV+NEDEV+NVREV),
!    INFOEV(1:4,1:NFCEV+NEDEV), the output from routine PRMDF3.
!
!    Input, integer ( kind = 4 ) IVRT(1:*), real ( kind = 8 ) FACVAL(1:*), EDGVAL(1:*),
!    VRTVAL(1:*), arrays output from routine DSMDF3.
!
!    Input, real ( kind = 8 ) VCL(1:3,1:*), the vertex coordinate list.
!
!    Input, integer ( kind = 4 ) FACEP(1:3,1:*), the face pointer list.
!
!    Input, integer ( kind = 4 ) FVL(1:6,1:*), the face vertex list.
!
!    Input, integer ( kind = 4 ) PFL(1:2,1:*), the list of signed face indices for each
!    polyhedron.
!
!    Input, real ( kind = 8 ) CNTR(1:3), the weighted centroid of polyhedron.
!
!    Input, integer ( kind = 4 ) NF, the number of faces in polyhedron if 1 < NF,
!    else NF = 1; in former case INDF, MEANF, STDVF values may be computed.
!
!    Output, real ( kind = 8 ) MDFINT, the integral of mdf in polyhedron.
!
!    Output, real ( kind = 8 ) MEAN, the mean mdf value in polyhedron.
!
!    Output, real ( kind = 8 ) STDV, the standard deviation of mdf in
!    polyhedron.
!
!    Output, real ( kind = 8 ) VOLP, the volume of polyhedron.
!
!    Output, integer ( kind = 4 ) INDF(1:NF), the indices in FACEP of faces of polyhedron.
!
!    Output, real ( kind = 8 ) MEANF(1:NF), the mean mdf value associated
!    with faces of polyhedron.
!
!    Output, real ( kind = 8 ) STDVF(1:NF), the standard deviation of
!    mdf associated with faces.
!
!    [Note: Above 3 arrays are computed if 1 < NF and (HFLAG =
!    FALSE or 0 < NFCEV + NEDEV + NVREV].
!
  implicit none

  integer ( kind = 4 ) nf
  integer ( kind = 4 ), parameter :: nqpt = 4

  real    ( kind = 8 ) cntr(3)
  real    ( kind = 8 ) cntrf(3)
  integer ( kind = 4 ), parameter :: edga = 5
  integer ( kind = 4 ), parameter :: edgc = 6
  real    ( kind = 8 ) edgval(*)
  logical              eval
  integer ( kind = 4 ) f
  integer ( kind = 4 ) facep(3,*)
  integer ( kind = 4 ), parameter :: facn = 2
  real    ( kind = 8 ) facval(*)
  integer ( kind = 4 ) fvl(6,*)
  logical              hflag
  integer ( kind = 4 ) headp
  integer ( kind = 4 ) i
  integer ( kind = 4 ) indf(nf)
  real    ( kind = 8 ) infoev(4,*)
  integer ( kind = 4 ) ivrt(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) lj
  integer ( kind = 4 ) lk
  integer ( kind = 4 ) listev(*)
  integer ( kind = 4 ), parameter :: loc = 1
  integer ( kind = 4 ) m
  real    ( kind = 8 ) mdf3
  real    ( kind = 8 ) mdfif
  real    ( kind = 8 ) mdfint
  real    ( kind = 8 ) mdfsf
  real    ( kind = 8 ) mdfsqi
  real    ( kind = 8 ) mean
  real    ( kind = 8 ) meanf(nf)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nedev
  integer ( kind = 4 ) nfcev
  integer ( kind = 4 ) nvrev
  integer ( kind = 4 ) pfl(2,*)
  integer ( kind = 4 ), parameter :: pred = 4
  integer ( kind = 4 ) q
  real    ( kind = 8 ), parameter, dimension ( 4, nqpt ) :: qc = reshape ( &
  (/ 0.5854102D+00, 0.1381966D+00, 0.1381966D+00, 0.1381966D+00, &
     0.1381966D+00, 0.5854102D+00, 0.1381966D+00, 0.1381966D+00, &
     0.1381966D+00, 0.1381966D+00, 0.5854102D+00, 0.1381966D+00, &
     0.1381966D+00, 0.1381966D+00, 0.1381966D+00, 0.5854102D+00/), &
    (/ 4, nqpt /) )
  real    ( kind = 8 ) stdv
  real    ( kind = 8 ) stdvf(nf)
  integer ( kind = 4 ), parameter :: succ = 3
  real    ( kind = 8 ) sum1
  real    ( kind = 8 ) sum2
  real    ( kind = 8 ) temp
  real    ( kind = 8 ) umdf
  real    ( kind = 8 ) val
  real    ( kind = 8 ) vcl(3,*)
  real    ( kind = 8 ) volf
  real    ( kind = 8 ) volp
  real    ( kind = 8 ) volt
  real    ( kind = 8 ) volth
  real    ( kind = 8 ) vrtval(*)
  real    ( kind = 8 ) widp
  real    ( kind = 8 ), parameter, dimension ( nqpt ) :: wt = (/ &
    0.25D+00, 0.25D+00, 0.25D+00, 0.25D+00 /)
  real    ( kind = 8 ) x
  real    ( kind = 8 ) y
  real    ( kind = 8 ) z
!
!  NQPT is number of quad points for numerical integration in tetrahedron
!  WT(I) is weight of Ith quadrature point
!  QC(1:4,I) are barycentric coordinates of Ith quadrature point
!
  if ( hflag ) then
    eval = ( 0 < nfcev + nedev + nvrev )
  else
    eval = .true.
  end if

  mdfint = 0.0D+00
  mdfsqi = 0.0D+00
  volp = 0.0D+00
  m = 0
  i = headp

  do

    f = abs ( pfl(1,i) )
    n = 0
    cntrf(1:3) = 0.0D+00
    j = facep(1,f)

    do

      lj = fvl(loc,j)
      cntrf(1:3) = cntrf(1:3) + vcl(1:3,lj)
      n = n + 1
      j = fvl(succ,j)

      if ( j == facep(1,f) ) then
        exit
      end if

    end do

    cntrf(1:3) = cntrf(1:3) / real ( n, kind = 8 )
    mdfif = 0.0D+00
    mdfsf = 0.0D+00
    volf = 0.0D+00
    lj = fvl(loc,j)

30  continue

    k = fvl(succ,j)
    lk = fvl(loc,k)
    volt = volth ( cntr, cntrf, vcl(1,lj), vcl(1,lk) ) / 6.0D+00
    volf = volf + volt

    if ( eval ) then

      sum1 = 0.0D+00
      sum2 = 0.0D+00

      do q = 1,nqpt

        x = qc(1,q)*cntr(1) + qc(2,q)*cntrf(1) + &
            qc(3,q)*vcl(1,lj) + qc(4,q)*vcl(1,lk)
        y = qc(1,q)*cntr(2) + qc(2,q)*cntrf(2) + &
            qc(3,q)*vcl(2,lj) + qc(4,q)*vcl(2,lk)
        z = qc(1,q)*cntr(3) + qc(2,q)*cntrf(3) + &
            qc(3,q)*vcl(3,lj) + qc(4,q)*vcl(3,lk)

        if ( hflag ) then
          val = mdf3(x,y,z,widp,nfcev,nedev,nvrev,listev, &
            infoev,ivrt,facval,edgval,vrtval,vcl)
        else
          val = umdf(x,y,z)
        end if

        temp = wt(q)*val
        sum1 = sum1 + temp
        sum2 = sum2 + temp*val

      end do

      mdfif = mdfif + sum1*volt
      mdfsf = mdfsf + sum2*volt

    end if

    j = k
    lj = lk

    if ( j /= facep(1,f)) go to 30

    if ( eval ) then

      mdfint = mdfint + mdfif
      mdfsqi = mdfsqi + mdfsf

      if ( 1 < nf ) then
        m = m + 1
        indf(m) = f
        meanf(m) = mdfif/volf
        temp = mdfsf/volf - meanf(m)**2
        stdvf(m) = sqrt ( max ( temp, 0.0D+00 ) )
      end if

    end if

    volp = volp + volf
    i = pfl(2,i)

    if ( i == headp ) then
      exit
    end if

  end do

  if ( eval ) then
    mean = mdfint / volp
    stdv = mdfsqi / volp - mean**2
    stdv = sqrt ( max ( stdv, 0.0D+00 ) )
  else
    mean = 1.0D+00 / widp**3
    mdfint = mean * volp
    stdv = 0.0D+00
  end if

  return
end
