!============================================================================
!                             April 10, 2002, v6.01
!                             February 23, 2003, v6.1

!   Ref[1]: "New Generation of Parton Distributions with Uncertainties from Global QCD Analysis"
!       By: J. Pumplin, D.R. Stump, J.Huston, H.L. Lai, P. Nadolsky, W.K. Tung
!       JHEP 0207:012(2002), hep-ph/0201195

!   Ref[2]: "Inclusive Jet Production, Parton Distributions, and the Search for New Physics"
!       By : D. Stump, J. Huston, J. Pumplin, W.K. Tung, H.L. Lai, S. Kuhlmann, J. Owens
!       hep-ph/0303013

!   This package contains
!   (1) 4 standard sets of CTEQ6 PDF's (CTEQ6M, CTEQ6D, CTEQ6L, CTEQ6L1) ;
!   (2) 40 up/down sets (with respect to CTEQ6M) for uncertainty studies from Ref[1];
!   (3) updated version of the above: CTEQ6.1M and its 40 up/down eigenvector sets from Ref[2].

!  The CTEQ6.1M set provides a global fit that is almost equivalent in every respect
!  to the published CTEQ6M, Ref[1], although some parton distributions (e.g., the gluon)
!  may deviate from CTEQ6M in some kinematic ranges by amounts that are well within the
!  specified uncertainties.
!  The more significant improvements of the new version are associated with some of the
!  40 eigenvector sets, which are made more symmetrical and reliable in (3), compared to (2).

!  Details about calling convention are:
! ---------------------------------------------------------------------------
!  Iset   PDF-set     Description       Alpha_s(Mz)**Lam4  Lam5   Table_File
! ===========================================================================
! Standard, "best-fit", sets:
! --------------------------
!   1    CTEQ6M   Standard MSbar scheme   0.118     326   226    cteq6m.tbl
!   2    CTEQ6D   Standard DIS scheme     0.118     326   226    cteq6d.tbl
!   3    CTEQ6L   Leading Order           0.118**   326** 226    cteq6l.tbl
!   4    CTEQ6L1  Leading Order           0.130**   215** 165    cteq6l1.tbl
! ============================================================================
! For uncertainty calculations using eigenvectors of the Hessian:
! ---------------------------------------------------------------
!     central + 40 up/down sets along 20 eigenvector directions
!                             -----------------------------
!                Original version, Ref[1]:  central fit: CTEQ6M (=CTEQ6M.00)
!                             -----------------------
!  1xx  CTEQ6M.xx  +/- sets               0.118     326   226    cteq6m1xx.tbl
!        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
!        e.g. 100      is CTEQ6M.00 (=CTEQ6M),
!             101/102 are CTEQ6M.01/02, +/- sets of 1st eigenvector, ... etc.
!                              -----------------------
!                Updated version, Ref[2]:  central fit: CTEQ6.1M (=CTEQ61.00)
!                              -----------------------
!  2xx  CTEQ61.xx  +/- sets               0.118     326   226    ctq61.xx.tbl
!        where xx = 01-40: 01/02 corresponds to +/- for the 1st eigenvector, ... etc.
!        e.g. 200      is CTEQ61.00 (=CTEQ6.1M),
!             201/202 are CTEQ61.01/02, +/- sets of 1st eigenvector, ... etc.
! ===========================================================================
!   ** ALL fits are obtained by using the same coupling strength
!   \alpha_s(Mz)=0.118 and the NLO running \alpha_s formula, except CTEQ6L1
!   which uses the LO running \alpha_s and its value determined from the fit.
!   For the LO fits, the evolution of the PDF and the hard cross sections are
!   calculated at LO.  More detailed discussions are given in the references.

!   The table grids are generated for 10^-6 < x < 1 and 1.3 < Q < 10,000 (GeV).
!   PDF values outside of the above range are returned using extrapolation.
!   Lam5 (Lam4) represents Lambda value (in MeV) for 5 (4) flavors.
!   The matching alpha_s between 4 and 5 flavors takes place at Q=4.5 GeV,
!   which is defined as the bottom quark mass, whenever it can be applied.

!   The Table_Files are assumed to be in the working directory.

!   Before using the PDF, it is necessary to do the initialization by
!       Call SetCtq6(Iset)
!   where Iset is the desired PDF specified in the above table.

!   The function Ctq6Pdf (Iparton, X, Q)
!   returns the parton distribution inside the proton for parton [Iparton]
!   at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset].
!   Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5)
!                            for (b, c, s, d, u, g, u_bar, ..., b_bar),

!   For detailed information on the parameters used, e.q. quark masses,
!   QCD Lambda, ... etc.,  see info lines at the beginning of the
!   Table_Files.

!   These programs, as provided, are in double precision.  By removing the
!   "Implicit Double Precision" lines, they can also be run in single
!   precision.

!   If you have detailed questions concerning these CTEQ6 distributions,
!   or if you find problems/bugs using this package, direct inquires to
!   Pumplin@pa.msu.edu or Tung@pa.msu.edu.

!===========================================================================

    Function Ctq6Pdf (Iparton, X, Q)
    Implicit Double Precision (A-H,O-Z)
    Logical :: Warn
    Common &
    / CtqPar2 / Nx, Nt, NfMx &
    / QCDtable /  Alambda, Nfl, Iorder

    Data Warn / .TRUE. /
    save Warn

    If (X < 0D0 .OR. X > 1D0) Then
        Print *, 'X out of range in Ctq6Pdf: ', X
    ! cccccccccccccccccc SM
        Ctq6Pdf = 0D0
        Return
    !        Stop
    ! ccccccccccccccccccccc
    Endif
    If (Q < Alambda) Then
        Print *, 'Q out of range in Ctq6Pdf: ', Q
    ! cccccccccccccccccc SM
        Ctq6Pdf = 0D0
        Return
    !        Stop
    ! ccccccccccccccccccccc
    Endif
    If ((Iparton < -NfMx .OR. Iparton > NfMx)) Then
        If (Warn) Then
        !        put a warning for calling extra flavor.
            Warn = .FALSE. 
            Print *, 'Warning: Iparton out of range in Ctq6Pdf! '
            Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx
        Endif
        Ctq6Pdf = 0D0
        Return
    Endif

    Ctq6Pdf = PartonX6 (Iparton, X, Q)
    if(Ctq6Pdf < 0.D0)  Ctq6Pdf = 0.D0

    Return

!                             ********************
    end Function Ctq6Pdf

    Subroutine SetCtq6 (Iset)
    Implicit Double Precision (A-H,O-Z)
    Parameter (Isetmax0=5)
    Character Flnm(Isetmax0)*6, nn*3, Tablefile*40
    Data (Flnm(I), I=1,Isetmax0) &
    / 'cteq6m', 'cteq6d', 'cteq6l', 'cteq6l','ctq61.'/
    Data Isetold, Isetmin0, Isetmin1, Isetmax1 /-987,1,100,140/
    Data Isetmin2,Isetmax2 /200,240/
    save

!             If data file not initialized, do so.
    If(Iset /= Isetold) then
        IU= NextUn()
        If (Iset >= Isetmin0 .AND. Iset <= 3) Then
            Tablefile=Flnm(Iset)//'.tbl'
        Elseif (Iset == 4) Then
            Tablefile=Flnm(Iset)//'1.tbl'
        Elseif (Iset >= Isetmin1 .AND. Iset <= Isetmax1) Then
            write(nn,'(I3)') Iset
            Tablefile=Flnm(1)//nn//'.tbl'
        Elseif (Iset >= Isetmin2 .AND. Iset <= Isetmax2) Then
            write(nn,'(I3)') Iset
            Tablefile=Flnm(5)//nn(2:3)//'.tbl'
        Else
            Print *, 'Invalid Iset number in SetCtq6 :', Iset
            Stop
        Endif
        Open(IU, File='pdfs/'//Tablefile, Status='OLD', Err=100)
        21 Call ReadTbl (IU)
        Close (IU)
        Isetold=Iset
    Endif
    Return

    100 Print *, ' Data file ', Tablefile, ' cannot be opened ' &
    //'in SetCtq6!!'
    Stop
!                             ********************
    end Subroutine SetCtq6

    Subroutine ReadTbl (Nu)
    Implicit Double Precision (A-H,O-Z)
    Character Line*80
    PARAMETER (MXX = 96, MXQ = 20, MXF = 5)
    PARAMETER (MXPQX = (MXF + 3) * MXQ * MXX)
    Common &
    / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX) &
    / CtqPar2 / Nx, Nt, NfMx &
    / XQrange / Qini, Qmax, Xmin &
    / QCDtable /  Alambda, Nfl, Iorder &
    / Masstbl / Amass(6)

    Read  (Nu, '(A)') Line
    Read  (Nu, '(A)') Line
    Read  (Nu, *) Dr, Fl, Al, (Amass(I),I=1,6)
    Iorder = Nint(Dr)
    Nfl = Nint(Fl)
    Alambda = Al

    Read  (Nu, '(A)') Line
    Read  (Nu, *) NX,  NT, NfMx

    Read  (Nu, '(A)') Line
    Read  (Nu, *) QINI, QMAX, (TV(I), I =0, NT)

    Read  (Nu, '(A)') Line
    Read  (Nu, *) XMIN, (XV(I), I =0, NX)

    Do 11 Iq = 0, NT
        TV(Iq) = Log(Log (TV(Iq) /Al))
    11 END DO

!                  Since quark = anti-quark for nfl>2 at this stage,
!                  we Read  out only the non-redundent data points
!     No of flavors = NfMx (sea) + 1 (gluon) + 2 (valence)

    Nblk = (NX+1) * (NT+1)
    Npts =  Nblk  * (NfMx+3)
    Read  (Nu, '(A)') Line
    Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts)

    Return
!                        ****************************
    end Subroutine ReadTbl

    Function NextUn()
!                                 Returns an unallocated FORTRAN i/o unit.
    Logical :: EX

    Do 10 N = 10, 300
        INQUIRE (UNIT=N, OPENED=EX)
        If ( .NOT. EX) then
            NextUn = N
            Return
        Endif
    10 END DO
    Stop ' There is no available I/O unit. '
!               *************************
    end Function NextUn


    SUBROUTINE POLINT (XA,YA,N,X,Y,DY)

    IMPLICIT DOUBLE PRECISION (A-H, O-Z)
!                                        Adapted from "Numerical Recipes"
    PARAMETER (NMAX=10)
    DIMENSION XA(N),YA(N),C(NMAX),D(NMAX)
    NS=1
    DIF=ABS(X-XA(1))
    DO 11 I=1,N
        DIFT=ABS(X-XA(I))
        IF (DIFT < DIF) THEN
            NS=I
            DIF=DIFT
        ENDIF
        C(I)=YA(I)
        D(I)=YA(I)
    11 END DO
    Y=YA(NS)
    NS=NS-1
    DO 13 M=1,N-1
        DO 12 I=1,N-M
            HO=XA(I)-X
            HP=XA(I+M)-X
            W=C(I+1)-D(I)
            DEN=HO-HP
        !          IF(DEN.EQ.0.)PAUSE
            DEN=W/DEN
            D(I)=HP*DEN
            C(I)=HO*DEN
        12 END DO
        IF (2*NS < N-M)THEN
            DY=C(NS+1)
        ELSE
            DY=D(NS)
            NS=NS-1
        ENDIF
        Y=Y+DY
    13 END DO
    RETURN
    END SUBROUTINE POLINT

    Function PartonX6 (IPRTN, XX, QQ)

!  Given the parton distribution function in the array U in
!  COMMON / PEVLDT / , this routine interpolates to find
!  the parton distribution at an arbitray point in x and q.

    Implicit Double Precision (A-H,O-Z)

    Parameter (MXX = 96, MXQ = 20, MXF = 5)
    Parameter (MXQX= MXQ * MXX,   MXPQX = MXQX * (MXF+3))

    Common &
    / CtqPar1 / Al, XV(0:MXX), TV(0:MXQ), UPD(MXPQX) &
    / CtqPar2 / Nx, Nt, NfMx &
    / XQrange / Qini, Qmax, Xmin

    Dimension fvec(4), fij(4)
    Dimension xvpow(0:mxx)
    Data OneP / 1.00001 /
    Data xpow / 0.3d0 /       !**** choice of interpolation variable
    Data nqvec / 4 /
    Data ientry / 0 /
    Save ientry,xvpow

! store the powers used for interpolation on first call...
    if(ientry == 0) then
        ientry = 1

        xvpow(0) = 0D0
        do i = 1, nx
            xvpow(i) = xv(i)**xpow
        enddo
    endif

    X = XX
    Q = QQ
    tt = log(log(Q/Al))

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
        Print '(A,1pE12.4)', 'Severe error: x <= 0 in PartonX6! x = ', x
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
        Print '(A,1pE12.4)', 'Severe error: x > 1 in PartonX6! x = ', x
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

    If (Iprtn >= 3) Then
        Ip = - Iprtn
    Else
        Ip = Iprtn
    EndIf
    jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1

    Do it = 1, nqvec

        J1  = jtmp + it*(NX+1)

        If (Jx == 0) Then
        !                          For the first 4 x points, interpolate x^2*f(x,Q)
        !                           This applies to the two lowest bins JLx = 0, 1
        !            We can not put the JLx.eq.1 bin into the "interrior" section
        !                           (as we do for q), since Upd(J1) is undefined.
            fij(1) = 0
            fij(2) = Upd(J1+1) * XV(1)**2
            fij(3) = Upd(J1+2) * XV(2)**2
            fij(4) = Upd(J1+3) * XV(3)**2
        
        !                 Use Polint which allows x to be anywhere w.r.t. the grid

            Call Polint (XVpow(0), Fij(1), 4, ss, Fx, Dfx)

            If (x > 0D0)  Fvec(it) =  Fx / x**2
        !                                              Pdf is undefined for x.eq.0
        ElseIf  (JLx == Nx-1) Then
        !                                                This is the highest x bin:

            Call Polint (XVpow(Nx-3), Upd(J1), 4, ss, Fx, Dfx)

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
        Call Polint (TV(0), Fvec(1), 4, tt, ff, Dfq)

    ElseIf (JLq >= Nt-1) Then
    !                         Last Q-bin, as well as extrapolation to higher Q
        Call Polint (TV(Nt-3), Fvec(1), 4, tt, ff, Dfq)
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

    PartonX6 = ff

    Return
!                                       ********************
    end Function PartonX6
