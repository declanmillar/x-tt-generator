!============================================================================
!                CTEQ-TEA Parton Distribution Functions: version 2014
!                       December 31, 2014

!   When using these PDFs, please cite the references below


!   This package provides a standard interface for CT10, CT12
!   (unpublished), and CT14 parton distribution functions.

!   The following sets of CTEQ PDF table files can be computed
!    with this program:
!               PDF                             References
!   (1) 1+50 sets of CT10 NNLO PDF's;             [1]
!   (2' ) 1+52 sets of CT10 NLO PDF's;            [2]
!   (2'') 1+52 sets of CT10W NLO PDF's;           [2]
!   (3) 4 sets of CT10W NNLO and NLO PDF's        [2]
!       with alternative alpha_s values;
!   (4) 1+56 sets of CT14 NNLO PDF's;             [3]
!   (5) 1+56 sets of CT14 NLO PDF's;              [3]
!   (6) 2 sets of CT14 LO PDF's;                  [3]
!   (7) 11 CT14 NNLO sets and 11 CT14 NLO sets    [3]
!       with alternative alpha_s values
!   (8) 3 CT14 NNLO and 3 CT14 NLO sets           [3]
!       with up to 3, 4, and 6 active quark flavors
!   (9) 4 CT14 NNLO sets with intrinsic charm     [X]

!   References
!   [1] J. Gao, M. Guzzi, J. Huston, H.-L. Lai, Z. Li, P. M. Nadolsky,
!       J. Pumplin, D. Stump, C.-P. Yuan,  arXiv:1302.6246 [hep-ph]
!   [2] H.-L. Lai, M. Guzzi, J. Huston, Z. Li, P. M. Nadolsky,
!       J. Pumplin, and C.-P. Yuan, arXiv: 1007.2241 [hep-ph]
!   [3] S. Dulat, T.-J. Hou, J. Gao, M. Guzzi, J. Huston,
!       P. M. Nadolsky, J. Pumplin, C. Schmidt, D. Stump, and
!       C.-P. Yuan, arXiv:1506.XXXX

! ===========================================================================
!   The table grids are generated for
!    *  10^-9 < x < 1 and 1.3 < Q < 10^5 (GeV).

!   PDF values outside of the above range are returned using extrapolation.

!   The Table_Files are assumed to be in the working directory.

!   Before using the PDF, it is necessary to do the initialization by
!       Call SetCT14(Iset)
!   where Tablefile is a 40-character text string with the name of the
!   the desired PDF specified in the above table.table (.pds) file

!   Other provided functions include:
!   The function CT14Pdf (Iparton, X, Q)
!     returns the parton distribution inside the proton for parton [Iparton]
!     at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
!     Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
!                              for (b, c, s, d, u, g, u_bar, ..., b_bar),

!   The function CT14Alphas (Q)
!     returns the value of the QCD coupling strength alpha_s(Q) at
!     an energy scale Q. The alpha_s value is obtained from the interpolation
!     of a numerical table of alpha_s included in each .pds file.
!     In the PDF fit, the alpha_s values are obtained at each momentum
!     scale by evolution in the HOPPET program at the respective QCD order
!     (NLO or NNLO). The table of alpha_s values at discrete Q values
!     is included in the input .pds file. The function CT14Alphas
!     estimates alpha_s at an arbitrary Q value, which agrees
!     with the direct evolution by HOPPET within a fraction of percent
!     point at typical Q.

!   The function CT14Mass(i)
!     returns the value of the quark mass for the i-th flavor.
!     The flavors are:
!     1  2  3  4  5  6
!     u  d  s  c  b  t

!   Values of various PDF parameters assumed in the computation of the
!    PDFs can be obtained by
!     Call CT14GetPars( xmin,Qini,Qmax,Nloops,Nfl),
!   which returns
!     xmin, the minimal value of x;
!     Qmin,  the initial Q scale for the PDF evolution;
!     Qmax,  the maximal Q scale included in the PDF table;
!     Nloop, the number of QCD loops (order of the PDF in the QCD coupling);
!     Nfl,   the maximal number of quark flavors assumed in the PDF and
!            alpha_s evolution.

!   These programs, as provided, are in double precision.  By removing the
!   "Implicit Double Precision" lines, they can also be run in single
!   precision.
!   If you have detailed questions concerning these CT14 distributions,
!   or if you find problems/bugs using this package, direct inquires to
!   nadolsky@smu.edu.

!===========================================================================

    Function CT14Pdf (Iparton, X, Q)
    Implicit Double Precision (A-H,O-Z)
    Logical :: Warn
    integer :: isetch, ipdsformat
    Common &
    / CtqPar2 / Nx, Nt, NfMx, MxVal &
    / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl  & !for external use
    /Setchange/ Isetch, ipdsset, ipdsformat

    Data Warn / .TRUE. /
    Data Qsml /.3d0/
    save Warn

    if (ipdsset /= 1) &
    STOP 'CT14Pdf: the PDF table was not initialized'

    If (X < 0d0 .OR. X > 1D0) Then
        Print *, 'X out of range in CT14Pdf: ', X
        CT14Pdf = 0D0
        Return
    Endif

    If (Q < Qsml) Then
        Print *, 'Q out of range in CT14Pdf: ', Q
        Stop
    Endif

    If (abs(Iparton) > NfMx) Then
        If (Warn) Then
        !        print a warning for calling extra flavor
            Warn = .FALSE.
            Print *, 'Warning: Iparton out of range in CT14Pdf! '
            Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
        Endif
        CT14Pdf = 0D0
    else

        CT14Pdf = PartonX12 (Iparton, X, Q)
        if (CT14Pdf < 0D0) CT14Pdf = 0D0
    endif                     !if (abs(Iparton...

    Return

!                             ********************
    end Function CT14Pdf

    Subroutine SetCT14(Tablefile)
    Implicit Double Precision (A-H,O-Z)
    Character Tablefile*40
    Common /Setchange/ Isetch, ipdsset, ipdsformat
    data ipdsset, ipdsformat/0,0/
    save

    IU= NextUn14()
    Open(IU, File='pdfs/'//Tablefile, Status='OLD', Err=100)
    Call Readpds0 (IU)
    Close (IU)
    Isetch=1; ipdsset=1
    Return

    100 Print *, ' Data file ', Tablefile, ' cannot be opened ' &
    //'in SetCT14!!'
    Stop
!                             ********************
    end Subroutine SetCT14

    subroutine CT14GetPars(xmin,Qini,Qmax,Nloops,Nfl)
! Get various parameters associated with the PDF grid
! Output: xmin  is the minimal value of x
!         Qmin  is the initial Q scale
!         Qmax  is the maximal Q scale
!         Nloop is the number of QCD loops
!         Nfl   is the maximal number of quark flavors
    implicit none
    double precision :: Qini0, Qmax0, Xmin0, xmin, Qini, Qmax
    integer :: Nloops, Ipk, Iorder, Nfl,Nfl0
    double precision :: AlfaQ, Qalfa

    common / XQrange / Qini0, Qmax0, Xmin0
    common / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl0

    Qini=Qini0; Qmax=Qmax0; Xmin=Xmin0
    Nloops=Iorder-1; Nfl=Nfl0

    return
    end subroutine CT14GetPars


    Function CT14Alphas (QQ)

    Implicit Double Precision (A-H,O-Z)

    PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
    PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
    double precision :: Alsout

    Common &
    / CtqPar1 / qBase,XV(0:MXX), TV(0:MXQ), UPD(MXPQX),AlsCTEQ(0:mxq) &
    / CtqPar2 / Nx, Nt, NfMx, MxVal &
    / XQrange / Qini, Qmax, Xmin &
    /Setchange/ Isetch, ipdsset, ipdsformat

    Data Q, JQ /-1D0, 0/
    save

    if (ipdsset /= 1) &
    STOP 'CT14Alphas: the PDF table was not initialized'


    if (ipdsformat < 11) then
        print *
        print *, &
        'STOP in CT14alphas: the PDF table file has an older format'
        print *, &
        'and does not include the table of QCD coupling values.'
        print *, &
        'You can still compute the PDFs, but do not call'
        print *, &
        'the CT14alphas function for the interpolation of alpha_s.'
        stop
    endif

    Q = QQ
    tt = log(log(Q/qBase))

!         --------------   Find lower end of interval containing Q, i.e.,
!                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
    JLq = -1
    JU = NT+1
    13 If (JU-JLq > 1) Then
        JM = (JU+JLq) / 2
        If (tt >= TV(JM)) Then
            JLq = JM
        Else
            JU = JM
        Endif
        Goto 13
    Endif

    If     (JLq <= 0) Then
        Jq = 0
    Elseif (JLq <= Nt-2) Then
    !                                  keep q in the middle, as shown above
        Jq = JLq - 1
    Else
    !                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

    Endif
!                                 This is the interpolation variable in Q
    Call Polint4F (TV(jq), AlsCTEQ(jq), tt, Alsout)

    CT14Alphas = Alsout

    Return
!                                       ********************
    end Function CT14Alphas


    function CT14Mass(i)
!     Returns the value of the quark mass for the i-th flavor
!     The flavors are:
!     1  2  3  4  5  6
!     u  d  s  c  b  t
    implicit none
    double precision :: CT14Mass, Amass
    integer ::  Isetch, ipdsset, i, ipdsformat
    common/Setchange/ Isetch, ipdsset, ipdsformat &
    / Masstbl / Amass(6)


    if (ipdsset /= 1) &
    STOP 'CT14Mass: the PDF table was not initialized'

    CT14Mass = Amass(i)

    return
    end function CT14Mass


    Subroutine Readpds0 (Nu)
    Implicit Double Precision (A-H,O-Z)
    Character Line*80
    integer :: ipdsformat
    PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
    PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)
    double precision :: qv(0:mxq)

    Common &
    / CtqPar1 / qBase,XV(0:MXX),TV(0:MXQ),UPD(MXPQX), AlsCTEQ(0:mxq) &
    / CtqPar2 / Nx, Nt, NfMx, MxVal &
    / XQrange / Qini, Qmax, Xmin &
    / Masstbl / Amass(6) &
    / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl  & !for external use
    /Setchange/ Isetch, ipdsset, ipdsformat

    Read  (Nu, '(A)') Line
    Read  (Nu, '(A)') Line

    if (Line(1:11) == '  ipk, Ordr') then !post-CT10 .pds format;
    ! Set alphas(MZ) at scale Zm, quark masses, and evolution type
        ipdsformat = 10           !Post-CT10 .pds format
        Read (Nu, *) ipk, Dr, Qalfa, AlfaQ, (amass(i),i=1,6)
        Iorder = Nint(Dr)
        read (Nu, '(A)') Line
        if (Line(1:7) == '  IMASS' ) then
            ipdsformat = 11         !CT12 .pds format
            read (Nu, *) aimass, fswitch, N0, N0, N0, Nfmx, MxVal
            Nfl=Nfmx
        else                      !Pre-CT12 format
            Read  (Nu, *) N0, N0, N0, NfMx, MxVal
        endif                     !Line(1:7)

    else                        !old .pds format;
        ipdsformat = 6            !CTEQ6.6 .pds format; alpha_s  is not specified
        Read (Nu, *) Dr, fl, Alambda, (amass(i),i=1,6)  !set Lambda_QCD
        Iorder = Nint(Dr); Nfl = Nint(fl)

        Read  (Nu, '(A)') Line
        Read  (Nu, *) dummy,dummy,dummy, NfMx, MxVal, N0
    endif                       !Line(1:11...

    Read  (Nu, '(A)') Line
    Read  (Nu, *) NX,  NT, N0, NG, N0

    if (ng > 0) Read  (Nu, '(A)') (Line, i=1,ng+1)

    Read  (Nu, '(A)') Line
    if (ipdsformat >= 11) then  !CT12 format with alpha_s values
        Read  (Nu, *) QINI, QMAX, (qv(i),TV(I), AlsCTEQ(I), I =0, NT)
    else                        !pre-CT12 format
        Read  (Nu, *) QINI, QMAX, (qv(i),TV(I), I =0, NT)
    endif                       !ipdsformat >= 11

! check that qBase is consistent with the definition of Tv(0:nQ) for 2 values of Qv
    qbase1 = Qv(1)/Exp(Exp(Tv(1)))
    qbase2 = Qv(nT)/Exp(Exp(Tv(NT)))
    if (abs(qbase1-qbase2) > 1e-5) then
        print *, 'Readpds0: something wrong with qbase'
        print *,'qbase1, qbase2=',qbase1,qbase2
        stop
    else
        qbase=(qbase1+qbase2)/2.0d0
    endif                     !abs(qbase1...

    Read  (Nu, '(A)') Line
    Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX)
    XV(0)=0D0

    Nblk = (NX+1) * (NT+1)
    Npts =  Nblk  * (NfMx+1+MxVal)
    Read  (Nu, '(A)') Line
    Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

    Return
!                        ****************************
    end Subroutine Readpds0

    Function PartonX12 (IPRTN, XX, QQ)

!  Given the parton distribution function in the array U in
!  COMMON / PEVLDT / , this routine interpolates to find
!  the parton distribution at an arbitray point in x and q.

    Implicit Double Precision (A-H,O-Z)

    PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4)
    PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX)

    Common &
    / CtqPar1 / qBase, XV(0:MXX), TV(0:MXQ),UPD(MXPQX),AlsCTEQ(0:mxq) &
    / CtqPar2 / Nx, Nt, NfMx, MxVal &
    / XQrange / Qini, Qmax, Xmin &
    /Setchange/ Isetch, ipdsset, ipdsformat

    Dimension fvec(4), fij(4)
    Dimension xvpow(0:mxx)
    Data OneP / 1.00001 /
    Data xpow / 0.3d0 /       !**** choice of interpolation variable
    Data nqvec / 4 /
    Data ientry / 0 /
    Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/
    Save xvpow
    Save X, Q, JX, JQ, JLX, JLQ
    Save ss, const1, const2, const3, const4, const5, const6
    Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3
    Save tmp1, tmp2, tdet

! store the powers used for interpolation on first call...
    if(Isetch == 1) then
        Isetch = 0

        xvpow(0) = 0D0
        do i = 1, nx
            xvpow(i) = xv(i)**xpow
        enddo
    elseIf((XX == X) .AND. (QQ == Q)) then
        goto 99
    endif

    X = XX
    Q = QQ
    tt = log(log(Q/qBase))

!      -------------    find lower end of interval containing x, i.e.,
!                       get jx such that xv(jx) .le. x .le. xv(jx+1)...
    JLx = -1
    JU = Nx+1
    11 If (JU-JLx > 1) Then
        JM = (JU+JLx) / 2
        If (X >= XV(JM)) Then
            JLx = JM
        Else
            JU = JM
        Endif
        Goto 11
    Endif
!                     Ix    0   1   2      Jx  JLx         Nx-2     Nx
!                           |---|---|---|...|---|-x-|---|...|---|---|
!                     x     0  Xmin               x                 1

    If     (JLx <= -1) Then
        Print '(A,1pE12.4)','Severe error: x <= 0 in PartonX12! x = ',x
        Stop
    ElseIf (JLx == 0) Then
        Jx = 0
    Elseif (JLx <= Nx-2) Then

    !                For interrior points, keep x in the middle, as shown above
        Jx = JLx - 1
    Elseif (JLx == Nx-1 .OR. x < OneP) Then

    !                  We tolerate a slight over-shoot of one (OneP=1.00001),
    !              perhaps due to roundoff or whatever, but not more than that.
    !                                      Keep at least 4 points >= Jx
        Jx = JLx - 2
    Else
        Print '(A,1pE12.4)','Severe error: x > 1 in PartonX12! x = ',x
        Stop
    Endif
!          ---------- Note: JLx uniquely identifies the x-bin; Jx does not.

!                       This is the variable to be interpolated in
    ss = x**xpow

    If (JLx >= 2 .AND. JLx <= Nx-2) Then

    !     initiation work for "interior bins": store the lattice points in s...
        svec1 = xvpow(jx)
        svec2 = xvpow(jx+1)
        svec3 = xvpow(jx+2)
        svec4 = xvpow(jx+3)

        s12 = svec1 - svec2
        s13 = svec1 - svec3
        s23 = svec2 - svec3
        s24 = svec2 - svec4
        s34 = svec3 - svec4

        sy2 = ss - svec2
        sy3 = ss - svec3

    ! constants needed for interpolating in s at fixed t lattice points...
        const1 = s13/s23
        const2 = s12/s23
        const3 = s34/s23
        const4 = s24/s23
        s1213 = s12 + s13
        s2434 = s24 + s34
        sdet = s12*s34 - s1213*s2434
        tmp = sy2*sy3/sdet
        const5 = (s34*sy2-s2434*sy3)*tmp/s12
        const6 = (s1213*sy2-s12*sy3)*tmp/s34

    EndIf

!         --------------Now find lower end of interval containing Q, i.e.,
!                          get jq such that qv(jq) .le. q .le. qv(jq+1)...
    JLq = -1
    JU = NT+1
    12 If (JU-JLq > 1) Then
        JM = (JU+JLq) / 2
        If (tt >= TV(JM)) Then
            JLq = JM
        Else
            JU = JM
        Endif
        Goto 12
    Endif

    If     (JLq <= 0) Then
        Jq = 0
    Elseif (JLq <= Nt-2) Then
    !                                  keep q in the middle, as shown above
        Jq = JLq - 1
    Else
    !                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq.
        Jq = Nt - 3

    Endif
!                                   This is the interpolation variable in Q

    If (JLq >= 1 .AND. JLq <= Nt-2) Then
    !                                        store the lattice points in t...
        tvec1 = Tv(jq)
        tvec2 = Tv(jq+1)
        tvec3 = Tv(jq+2)
        tvec4 = Tv(jq+3)

        t12 = tvec1 - tvec2
        t13 = tvec1 - tvec3
        t23 = tvec2 - tvec3
        t24 = tvec2 - tvec4
        t34 = tvec3 - tvec4

        ty2 = tt - tvec2
        ty3 = tt - tvec3

        tmp1 = t12 + t13
        tmp2 = t24 + t34

        tdet = t12*t34 - tmp1*tmp2

    EndIf


! get the pdf function values at the lattice points...

    99 If (Iprtn > MxVal) Then
        Ip = - Iprtn
    Else
        Ip = Iprtn
    EndIf
    jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

    Do it = 1, nqvec
        J1  = jtmp + it*(NX+1)

        If (Jx == 0) Then
        !                      For the first 4 x points, interpolate x^2*f(x,Q)
        !                      This applies to the two lowest bins JLx = 0, 1
        !            We can not put the JLx.eq.1 bin into the "interrior" section
        !                           (as we do for q), since Upd(J1) is undefined.
            fij(1) = 0
            fij(2) = Upd(J1+1) * XV(1)**2
            fij(3) = Upd(J1+2) * XV(2)**2
            fij(4) = Upd(J1+3) * XV(3)**2

        !                 Use Polint which allows x to be anywhere w.r.t. the grid

            Call Polint4F (XVpow(0), Fij(1), ss, Fx)

            If (x > 0D0)  Fvec(it) =  Fx / x**2
        !                                              Pdf is undefined for x.eq.0
        ElseIf  (JLx == Nx-1) Then
        !                                                This is the highest x bin:

            Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx)

            Fvec(it) = Fx

        Else
        !                       for all interior points, use Jon's in-line function
        !                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2)
            sf2 = Upd(J1+1)
            sf3 = Upd(J1+2)

            g1 =  sf2*const1 - sf3*const2
            g4 = -sf2*const3 + sf3*const4

            Fvec(it) = (const5*(Upd(J1)-g1) &
            + const6*(Upd(J1+3)-g4) &
            + sf2*sy3 - sf3*sy2) / s23

        Endif

    enddo
!                                   We now have the four values Fvec(1:4)
!     interpolate in t...

    If (JLq <= 0) Then
    !                         1st Q-bin, as well as extrapolation to lower Q
        Call Polint4F (TV(0), Fvec(1), tt, ff)

    ElseIf (JLq >= Nt-1) Then
    !                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint4F (TV(Nt-3), Fvec(1), tt, ff)
    Else
    !                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2)
    !       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for
    !                         the full range QV(0:Nt)  (in contrast to XV)
        tf2 = fvec(2)
        tf3 = fvec(3)

        g1 = ( tf2*t13 - tf3*t12) / t23
        g4 = (-tf2*t34 + tf3*t24) / t23

        h00 = ((t34*ty2-tmp2*ty3)*(fvec(1)-g1)/t12 &
        +  (tmp1*ty2-t12*ty3)*(fvec(4)-g4)/t34)

        ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23
    EndIf

    PartonX12 = ff

    Return
!                                       ********************
    end Function PartonX12


    SUBROUTINE POLINT4F (XA,YA,X,Y)

    IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes",
!  but assuming N=4, and ignoring the error estimation.
!  suggested by Z. Sullivan.
    DIMENSION XA(*),YA(*)

    H1=XA(1)-X
    H2=XA(2)-X
    H3=XA(3)-X
    H4=XA(4)-X

    W=YA(2)-YA(1)
    DEN=W/(H1-H2)
    D1=H2*DEN
    C1=H1*DEN

    W=YA(3)-YA(2)
    DEN=W/(H2-H3)
    D2=H3*DEN
    C2=H2*DEN

    W=YA(4)-YA(3)
    DEN=W/(H3-H4)
    D3=H4*DEN
    C3=H3*DEN

    W=C2-D1
    DEN=W/(H1-H3)
    CD1=H3*DEN
    CC1=H1*DEN

    W=C3-D2
    DEN=W/(H2-H4)
    CD2=H4*DEN
    CC2=H2*DEN

    W=CC2-CD1
    DEN=W/(H1-H4)
    DD1=H4*DEN
    DC1=H1*DEN

    If((H3+H4) < 0D0) Then
        Y=YA(4)+D3+CD2+DD1
    Elseif((H2+H3) < 0D0) Then
        Y=YA(3)+D2+CD1+DC1
    Elseif((H1+H2) < 0D0) Then
        Y=YA(2)+C2+CD1+DC1
    ELSE
        Y=YA(1)+C1+CC1+DC1
    ENDIF

    RETURN
!               *************************
    END SUBROUTINE POLINT4F

    Function NextUn14()
!                                 Returns an unallocated FORTRAN i/o unit.
    Logical :: EX

    Do 10 N = 10, 300
        INQUIRE (UNIT=N, OPENED=EX)
        If ( .NOT. EX) then
            NextUn14 = N
            Return
        Endif
    10 END DO
    Stop ' There is no available I/O unit. '
!               *************************
    end Function NextUn14
