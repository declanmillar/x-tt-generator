module lhef

    use vamp_kinds

    implicit none

    public :: lhe_open
    public :: lhe_beam
    public :: lhe_process
    public :: lhe_add_event
    public :: lhe_add_particle
    public :: lhe_close

    private :: mass

    integer, parameter, private :: lhe = 20
    integer, parameter, private :: npr = 1 ! the number of different user sub-processes (I fix at 1)

contains


function mass( p )

    real( kind = default ) :: mass, p(4)

    mass = sqrt( abs( p(4) * p(4) - p(1) * p(1) - p(2) * p(2) - p(3) * p(3) ) )

end function mass


subroutine lhe_open( lhe_file )

    character( * ) lhe_file

    open( unit = lhe, file = trim( lhe_file ), status = "replace", action = "write" )
    write( lhe, "( a )" ) '<LesHouchesEvents version="1.0">'
    write( lhe, "( a )" ) '<!--'
    write( lhe, "( a )" ) 'File generated with zprime-top-generator'
    write( lhe, "( a )" ) '-->'

end subroutine lhe_open


subroutine lhe_beam( idbm1, idbm2, ebm1, ebm2, pdfg1, pdfg2, pdfs1, pdfs2, idwt )

    integer :: idbm1, idbm2 ! ID of beam particle 1 and 2 according to the PDG
    real( kind = default ) :: ebm1, ebm2 ! energy in GeV of beam particles 1 and 2
    integer :: pdfg1, pdfg2 ! author group for beam 1 and 2 according to Cernlib PDFlib
    integer :: pdfs1, pdfs2 ! PDF set ID for beam 1 and 2 according to Cernlib PDFlib
    integer :: idwt ! master switch dictating how the event weights (XWGTUP) are interpreted

    write( lhe, "( a )" ) "<init>"
    write( lhe, "( i8, i8, eS14.6, eS14.6, i5, i5, i6, i6, i5, i5 )" ) &
        idbm1, idbm2, ebm1, ebm2, pdfg1, pdfg1, pdfs1, pdfs2, idwt, npr

end subroutine lhe_beam


subroutine lhe_process( xsec, xerr, xmax, lpr )

    real( kind = default ) :: xsec ! the cross section for process j in pb this entry is mandatory for idwt=±2.
    real( kind = default ) :: xerr ! the statistical error associated with the cross section of process j in pb
    real( kind = default ) :: xmax ! the maximum xwgt for process j
    integer :: lpr ! a list of all user process IDs that can appear in idpr of hepe for this run
    write( lhe, "( ES14.6, ES14.6, ES14.6, I5 )" ) xsec, xerr, xmax, lpr
    write( lhe, "( a )" ) "</init>"

end subroutine lhe_process


subroutine lhe_add_event( n, idpr, xwgt, scal, aqed, aqcd )

    integer :: n ! number of particle entries in this event
    integer :: idpr ! ID of the process for this event
    real( kind = default ) :: xwgt ! event weight
    real( kind = default ) :: scal ! scale of the event in GeV, as used for calculation of PDFs
    real( kind = default ) :: aqed ! the QED coupling alpha_QED used for this event (e.g. 1/128)
    real( kind = default ) :: aqcd ! the QCD coupling alpha_QCD used for this event

    write( lhe, "( a )" ) "<event>"
    write( lhe, "( i6, i6, es14.6, es14.6, es14.6, es14.6 )" ) n, idpr, xwgt, scal, aqed, aqcd

end subroutine lhe_add_event


subroutine lhe_add_particle( id, ist, moth1, moth2, icol1, icol2, p )

    integer :: id ! particle ID according to Particle Data Group convention
    integer :: ist ! status code
    integer :: moth1, moth2 ! index of first and last mother
    integer :: icol1, icol2 ! 1(2) is the integer tag for the (anti) color flow line passing through the color of the particle
    real( kind = default ) :: p(4) ! lab frame momentum
    real( kind = default ) :: vtim ! invariant lifetime ctau (distance from production to decay) in mm
    real( kind = default ) :: spin ! cosine of the angle between the spin-vector of particle I
                                   ! and the 3-momentum of the decaying particle, specified in the lab frame
    real( kind = default ) :: m

    vtim = 0.0
    spin = 9.0

    m = mass( p )
    if ( m < 1.d0 ) m = 0.d0

    write( lhe, "( i8, i5, i5, i5, i5, i5, es18.10, es18.10, es18.10, es18.10, es18.10, f3.0, f3.0)" ) &
        id, ist, moth1, moth2, icol1, icol2, p, m, vtim, spin

end subroutine lhe_add_particle


subroutine lhe_end_event

    write( lhe, "( a )" ) "</event>"

end subroutine lhe_end_event


subroutine lhe_close

    write( lhe, "( a )" ) "</LesHouchesEvents>"
    close( lhe )

end subroutine lhe_close


end module lhef
