module configuration

    use kinds

    implicit none

    ! provide config
    logical :: new_grid
    real(kind=default) :: sqrts
    integer :: ppbar
    integer :: final_state
    integer :: pdf
    character(50) :: model_name
    character(100) :: filename
    character(100) :: log_file
    character(100) :: lhe_file
    character(100) :: grid_file
    character(100) :: xsec_file
    logical :: include_signal
    logical :: include_background
    logical :: include_gg
    logical :: include_qq
    logical :: include_uu
    logical :: include_dd
    logical :: include_a
    logical :: include_z
    logical :: include_x
    integer :: interference
    logical :: nwa
    logical :: multichannel
    logical :: use_rambo
    logical :: flatten_integrand
    logical :: verbose
    real(kind=default) :: ecm_low, ecm_up
    integer :: ncall
    integer :: nevents
    integer :: itmx
    logical :: batch
    logical :: cut
    logical :: unweighted

    ! derived config
    integer :: nfinal
    integer :: tops_decay
    integer :: ffinal
    integer :: nloops
    integer :: z_mixing
    real(kind=default) :: s2mix
    integer :: pid
    integer :: idw
    integer :: pdf_group
    integer :: pdf_set

    ! constants
    real(kind=default), parameter :: pi = 3.14159265358979323846d0
    real(kind=default), parameter :: twopi = 6.283185307179586d0

    ! methods
    public :: read_config
    public :: modify_config
    public :: print_config

contains

subroutine read_config
    if (verbose) print*, "config: reading ..."
    read(5,*) new_grid
    read(5,"(a)") grid_file
    read(5,"(a)") xsec_file
    read(5,"(a)") log_file
    read(5,"(a)") lhe_file
    read(5,*) ppbar ! 0 = pp, 1 = ppbar
    read(5,*) final_state ! 1 = no decay, 1 = dilepton
    read(5,*) model_name
    read(5,*) pdf
    read(5,*) include_signal
    read(5,*) include_background
    read(5,*) include_gg
    read(5,*) include_qq
    read(5,*) include_uu
    read(5,*) include_dd
    read(5,*) include_a
    read(5,*) include_z
    read(5,*) include_x
    read(5,*) interference
    read(5,*) nwa ! 0 = actual top widths, 1 = tops in NWA
    read(5,*) sqrts
    read(5,*) itmx
    read(5,*) ncall
    read(5,*) nevents
    read(5,*) unweighted
    read(5,*) use_rambo
    read(5,*) flatten_integrand
    read(5,*) multichannel
    read(5,*) verbose
    read(5,*) ecm_low
    read(5,*) ecm_up
    read(5,*) batch
    read(5,*) cut
end subroutine read_config

subroutine modify_config

    integer :: i

    if (verbose) print*, "config: interpreting ..."

    if (final_state < 1) then
        nfinal = 4
    else
        nfinal = 8
    end if

    if (final_state == -1) then
        ffinal = 1
    else
        ffinal = 11
    end if

    if (final_state < 1) then
        tops_decay = 0
    else
        tops_decay = 1
    end if

    ! from http://lhapdf.hepforge.org/pdfsets.html
    ! https://arxiv.org/pdf/hep-ph/9907231.pdf
    ! www.hep.phy.cam.ac.uk/theory/webber/MCatNLO/pdflib_doc.ps.gz

    if      (pdf ==  1) then ! CTEQ6M
        pdf_group = 0
        pdf_set = 10000
    else if (pdf ==  2) then ! CTEQ6D
        pdf_group = 0
        pdf_set = 10000
    else if (pdf ==  3) then ! CTEQ6L
        pdf_group = 0
        pdf_set = 10000
    else if (pdf ==  4) then ! CTEQ6L1
        pdf_group = 0
        pdf_set = 10042
    else if (pdf ==  5) then ! MRS99 (cor01)
        pdf_group = 3
        pdf_set = 78
    else if (pdf ==  6) then ! MRS99 (cor02)
        pdf_group = 3
        pdf_set = 79
    else if (pdf ==  7) then ! MRS99 (cor03)
        pdf_group = 3
        pdf_set = 80
    else if (pdf ==  8) then ! MRS99 (cor04)
        pdf_group = 3
        pdf_set = 81
    else if (pdf ==  9) then ! MRS99 (cor05)
        pdf_group = 3
        pdf_set = 82
    else if (pdf == 10) then ! CT14LN
        pdf_group = 0
        pdf_set = 13200
    else if (pdf == 11) then ! CT14LL
        pdf_group = 0
        pdf_set = 13205
    end if

    if (unweighted) then
        idw = 3
    else
        idw = 2
    end if

end subroutine modify_config

subroutine print_config

    if (verbose) print*, "config: printing ..."

    if (verbose) print*, "output =", trim(lhe_file)
    if (verbose .and. new_grid) then
    print*, "output =", trim(grid_file)
    else if (verbose) then
        print*, "input = ", trim(grid_file)
    end if

    if (ppbar == 0) then
        if (final_state == -1) print*, "process = ", "p p -> l+ l-"
        if (final_state ==  0) print*, "process = ", "p p -> t t~"
        if (final_state ==  1) print*, "process = ", "p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~"
    else if (ppbar == 1) then
        if (final_state == -1) print*, "process = ", "p p~ -> l+ l-"
        if (final_state ==  0) print*, "process = ", "p p~ -> t t~"
        if (final_state ==  1) print*, "process = ", "p p~ -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~"
    end if
    print*, "model = ", model_name
    if (pdf ==  1) print*, "PDF = ", "CTEQ6M"
    if (pdf ==  2) print*, "PDF = ", "CTEQ6D"
    if (pdf ==  3) print*, "PDF = ", "CTEQ6L"
    if (pdf ==  4) print*, "PDF = ", "CTEQ6L1"
    if (pdf ==  5) print*, "PDF = ", "MRST 99 (g up)"
    if (pdf ==  6) print*, "PDF = ", "MRST 99 (g down)"
    if (pdf ==  7) print*, "PDF = ", "MRST 99 (g up)"
    if (pdf ==  8) print*, "PDF = ", "MRST 99 (g up)"
    if (pdf ==  9) print*, "PDF = ", "MRST 99 (g up)"
    if (pdf == 10) print*, "PDF = ", "CT14LN"
    if (pdf == 11) print*, "PDF = ", "CT14LL"
    print*, "unweighted = ", unweighted
    print*, "detector cuts = ", cut
    print*, "multichannel = ", multichannel
    print*, "flatten integrand = ", flatten_integrand
    if (ecm_low .ne. 0) print*, "Ecm low = ", ecm_low
    if (ecm_up .ne. 0) print*, "Ecm up = ", ecm_up
    print*, "NWA = ", nwa
    print*, "RAMBO = ", use_rambo
    print*, "include gg = ", include_gg
    print*, "include qq = ", include_qq
    print*, "include uu = ", include_uu
    print*, "include dd = ", include_dd
    print*, "include A = ", include_a
    print*, "include Z = ", include_z
    print*, "include Z' = ", include_x
    if (include_dd .and. include_uu) then
    if (interference == 0) print*, "interference = ", "(gamma) + (Z) + (Z')"
    if (interference == 1) print*, "interference = ", "(gamma + Z + Z')"
    if (interference == 2) print*, "interference = ", "(gamma + Z) + (Z')"
    if (interference == 3) print*, "interference = ", "(gamma + Z + Z') - (gamma) - (Z)"
    if (interference == 4) print*, "interference = ", "(gamma + Z + Z') - (gamma) - (Z) - (Z')"
    end if
    print*, "collider energy = ", sqrts
end subroutine print_config

end module configuration
