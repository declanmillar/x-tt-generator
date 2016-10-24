module configuration

  use integration, only: seed, itmx, ncall, acc

  implicit none

  ! config
  logical :: ntuple_out
  logical :: lhef_out
  real :: sqrts
  integer :: initial_state
  integer :: final_state
  integer :: ipdf
  character(50) :: model_name
  character(100) :: ntuple_file
  character(100) :: log_file
  character(100) :: lhe_file
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
  logical :: symmetrise

  logical :: use_rambo
  logical :: map_phase_space
  logical :: verbose
  real :: ecm_low, ecm_up
  logical :: cut

  ! derived config
  integer :: n_final
  integer :: tops_decay
  integer :: ffinal
  integer :: nloops
  real :: lambdaqcd4
  integer :: z_mixing
  real :: s2mix
  integer :: pid

  ! constants
  real, parameter :: pi = 3.14159265358979323846d0
  integer, parameter :: log = 10

  ! methods
  public :: read_config
  public :: modify_config
  public :: print_config

contains

  subroutine read_config

    ! read config file
    read(5,*) ntuple_out
    read(5,*) lhef_out
    read(5,"(a)") ntuple_file
    read(5,"(a)") lhe_file
    read(5,"(a)") log_file
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
    read(5,*) use_nwa ! 0:actual top widths,1: tops in nwa
    read(5,*) sqrts
    read(5,*) seed
    read(5,*) itmx
    read(5,*) ncall
    read(5,*) acc
    read(5,*) use_rambo
    read(5,*) map_phase_space
    read(5,*) symmetrise
    read(5,*) verbose
    read(5,*) ecm_low
    read(5,*) ecm_up
    read(5,*) cut

  end subroutine read_config

  subroutine modify_config

    integer :: i

    if (final_state <= 0) then
      n_final = 4
    else
      n_final = 8
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

  end subroutine modify_config

subroutine print_config

  if (initial_state == 0) then
    if (final_state == -1) then
      write(log,*) "Process:p p -> l+ l-"
    else if (final_state == 0) then
      write(log,*) 'Process:p p -> t t~'
    else if (final_state == 1) then
      write(log,*) 'Process:p p -> t t~ -> b b W+ W- -> b b l+ l- nu nu~'
    end if
  else if (initial_state == 1) then
    if (final_state == -1) then
      write(log,*) 'pp~ -> l+l-'
    else if (final_state == 0) then
      write(log,*) 'pp~ -> tt~'
    else if (final_state == 1) then
      write(log,*) 'Process:p p~ -> t t~ -> b b W+ W- -> b b l+ l- nu nu~'
    end if
  end if
  if (ipdf == 1) write(log,*) 'PDFs: CTEQ6m'
  if (ipdf == 2) write(log,*) 'PDFs: CTEQ6d'
  if (ipdf == 3) write(log,*) 'PDFs: CTEQ6l'
  if (ipdf == 4) write(log,*) 'PDFs: CTEQ6l1'
  if (ipdf == 5) write(log,*) 'PDFs: MRS99 (cor01)'
  if (ipdf == 6) write(log,*) 'PDFs: MRS99 (cor02)'
  if (ipdf == 7) write(log,*) 'PDFs: MRS99 (cor03)'
  if (ipdf == 8) write(log,*) 'PDFs: MRS99 (cor04)'
  if (ipdf == 9) write(log,*) 'PDFs: MRS99 (cor05)'
  if (ipdf == 10) write(log,*) 'PDFs: CT14LN'
  if (ipdf == 11) write(log,*) 'PDFs: CT14LL'
  if (final_state >= 1 .and. use_nwa) write(log,*) 'NWA:ON'
  if (final_state >= 1 .and. .not. use_nwa) write(log,*) 'NWA:OFF'
  write(log,*) 'Model:', model_name
  if (include_a) write(log,*) 'photon:ON '
  if (.not. include_a) write(log,*) 'photon:OFF'
  if (include_z) write(log,*) 'Z boson:ON '
  if (.not. include_z) write(log,*) 'Z boson:OFF'
  if (include_x) write(log,*) 'BSM:ON '
  if (.not. include_x) write(log,*) 'BSM:OFF'
  write(log,*) "include_gg: ", include_gg
  write(log,*) "include_qq: ", include_qq
  write(log,*) "include_uu: ", include_uu
  write(log,*) "include_dd: ", include_dd
  write(log,*) "include_a: ", include_a
  write(log,*) "include_z: ", include_z
  write(log,*) "include_x: ", include_x
  if (interference == 0) write(log,*) "Interference: none"
  if (interference == 1) write(log,*) "Interference: (gamma + Z) + (Z')"
  if (interference == 2) write(log,*) "Interference: (gamma + Z + Z')"
  if (interference == 3) write(log,*) "Interference: (gamma + Z + Z') - (gamma) - (Z)"
  if (interference == 4) write(log,*) "Interference: (gamma + Z + Z') - (gamma) - (Z) - (Z')"
  if (symmetrise) write(log,*) 'Symmetrising integration: x1<->x2!'
  if (use_rambo) write(log,*) 'RAMBO: ON'
  write(log,*) "map_phase_space: ", map_phase_space
  write(log,*) 'seed:', seed
  write(log,*) 'collider energy:', sqrts
  if (ecm_low > 0) write(log,*) "E_cm low         :", ecm_low
  if (ecm_up > 0) write(log,*) "E_cm up             ", ecm_up
  if (lhef_out) write(log,*) "Events written to:", lhe_file
  if (ntuple_out) write(log,*) "Events written to:", ntuple_file
end subroutine print_config

end module configuration
