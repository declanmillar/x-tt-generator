module configuration

  use kinds

  implicit none

  ! provide config
  logical :: ntuple_out
  logical :: lhef_out
  logical :: new_grid
  real(kind=default) :: sqrts
  integer :: initial_state
  integer :: final_state
  integer :: ipdf
  character(50) :: model_name
  character(100) :: ntuple_file
  character(100) :: log_file
  character(100) :: lhe_file
  character(100) :: grid_file
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
  logical :: use_nwa
  logical :: multichannel
  logical :: symmetrise
  logical :: use_rambo
  logical :: map_phase_space
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

  ! methods
  public :: read_config
  public :: modify_config
  public :: print_config

contains

subroutine read_config
  print*, "config: reading ..."
  read(5,*) ntuple_out
  read(5,*) lhef_out
  read(5,*) new_grid
  read(5,"(a)") ntuple_file
  read(5,"(a)") lhe_file
  read(5,"(a)") log_file
  read(5,"(a)") grid_file
  read(5,*) initial_state ! 0 = pp, 1 = ppbar
  read(5,*) final_state ! 1 = no decay, 1 = dilepton, 2 = semilepton, 4 = full hadron
  read(5,*) model_name
  read(5,*) ipdf
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
  read(5,*) use_nwa ! 0 = actual top widths, 1 = tops in NWA
  read(5,*) sqrts
  read(5,*) itmx
  read(5,*) ncall
  read(5,*) nevents
  read(5,*) unweighted
  read(5,*) use_rambo
  read(5,*) map_phase_space
  read(5,*) multichannel
  read(5,*) symmetrise
  read(5,*) verbose
  read(5,*) ecm_low
  read(5,*) ecm_up
  read(5,*) batch
  read(5,*) cut
end subroutine read_config

subroutine modify_config

  integer :: i

  print*, "config: interpreting ..."

  if (final_state <= 0) then
    nfinal = 4
  else
    nfinal = 8
  end if

  if (final_state == -1) then
    ffinal = 1
  else
    ffinal = 11
  end if

  if (final_state <= 0) then
    tops_decay = 0
  else
    tops_decay = 1
  end if

  ! from http://lhapdf.hepforge.org/pdfsets.html
  ! https://arxiv.org/pdf/hep-ph/9907231.pdf
  ! www.hep.phy.cam.ac.uk/theory/webber/MCatNLO/pdflib_doc.ps.gz

  if      (ipdf ==  1) then ! PDF: CTEQ6M
    pdf_group = 0
    pdf_set = 10000
  else if (ipdf ==  2) then ! PDF: CTEQ6D
    pdf_group = 0
    pdf_set = 10000
  else if (ipdf ==  3) then ! PDF: CTEQ6L
    pdf_group = 0
    pdf_set = 10000
  else if (ipdf ==  4) then ! PDF: CTEQ6L1
    pdf_group = 0
    pdf_set = 10042
  else if (ipdf ==  5) then ! PDF: MRS99 (cor01)
    pdf_group = 3
    pdf_set = 78
  else if (ipdf ==  6) then ! PDF: MRS99 (cor02)
    pdf_group = 3
    pdf_set = 79
  else if (ipdf ==  7) then ! PDF: MRS99 (cor03)
    pdf_group = 3
    pdf_set = 80
  else if (ipdf ==  8) then ! PDF: MRS99 (cor04)
    pdf_group = 3
    pdf_set = 81
  else if (ipdf ==  9) then ! PDF: MRS99 (cor05)
    pdf_group = 3
    pdf_set = 82
  else if (ipdf == 10) then ! PDF: CT14LN
    pdf_group = 0 
    pdf_set = 13200
  else if (ipdf == 11) then ! PDF: CT14LL
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

  print*, "config: printing ..."

  if (lhef_out) print*, "output: ", lhe_file
  if (ntuple_out) print*, "output: ", ntuple_file
  if (new_grid) then 
    print*, "output: ", grid_file
  else
    print*, "input:  ", grid_file
  end if

  if (initial_state == 0) then
    if (final_state == -1) then
      print*, "process: p p -> l+ l-"
    else if (final_state == 0) then
      print*, "process: p p -> t t~"
    else if (final_state == 1) then
      print*, "process: p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~"
    end if
  else if (initial_state == 1) then
    if (final_state == -1) then
      print*, "process: p p~ -> l+ l-"
    else if (final_state == 0) then
      print*, "process: p p~ -> t t~"
    else if (final_state == 1) then
      print*, "process: p p~ -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~"
    end if
  end if
  print*, "VAMP iterations: ", itmx
  print*, "VAMP calls: ", ncall
  print*, "number of events: ", nevents
  print*, "unweighted events: ", unweighted
  print*, "apply detector cuts: ", cut
  print*, "multichannel: ", multichannel
  print*, "symmetrise parton momentum fraction: ", symmetrise
  print*, "map phase space: ", map_phase_space
  print*, "collider energy: ", sqrts
  if (Ecm_low .ne. 0) print*, "Ecm low: ", ecm_low
  if (Ecm_up .ne. 0) print*, "Ecm up:  ", ecm_up
  print*, "NWA: ", use_nwa
  print*, "RAMBO: ", use_rambo
  print*, "model: ", model_name
  print*, "include gg: ", include_gg
  print*, "include qq: ", include_qq
  print*, "include uu: ", include_uu
  print*, "include dd: ", include_dd
  print*, "include A:  ", include_a
  print*, "include Z:  ", include_z
  print*, "include Z': ", include_x
  if (interference == 0) print*, "interference: (gamma) + (Z) + (Z')"
  if (interference == 1) print*, "interference: (gamma + Z + Z')"
  if (interference == 2) print*, "interference: (gamma + Z) + (Z')"
  if (interference == 3) print*, "interference: (gamma + Z + Z') - (gamma) - (Z)"
  if (interference == 4) print*, "interference: (gamma + Z + Z') - (gamma) - (Z) - (Z')"
  if (ipdf ==  1) print*, "PDF: CTEQ6M"
  if (ipdf ==  2) print*, "PDF: CTEQ6D"
  if (ipdf ==  3) print*, "PDF: CTEQ6L"
  if (ipdf ==  4) print*, "PDF: CTEQ6L1"
  if (ipdf ==  5) print*, "PDF: MRST 99 (g up)"
  if (ipdf ==  6) print*, "PDF: MRST 99 (g down)"
  if (ipdf ==  7) print*, "PDF: MRST 99 (g up)"
  if (ipdf ==  8) print*, "PDF: MRST 99 (g up)"
  if (ipdf ==  9) print*, "PDF: MRST 99 (g up)"
  if (ipdf == 10) print*, "PDF: CT14LN"
  if (ipdf == 11) print*, "PDF: CT14LL"
end subroutine print_config

end module configuration
