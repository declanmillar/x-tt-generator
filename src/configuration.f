module configuration

  use kinds

  implicit none

  ! provide config
  logical :: ntuple_out
  logical :: lhef_out
  real(kind=default) :: sqrts
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
  real(kind=default) :: ecm_low, ecm_up
  integer :: seed
  integer :: ncall
  integer :: nevents
  integer :: itmx

  ! derived config
  integer :: n_final
  integer :: tops_decay
  integer :: ffinal
  integer :: nloops
  real(kind=default) :: lambdaqcd4
  integer :: z_mixing
  real(kind=default) :: s2mix
  integer :: pid

  ! constants
  real(kind=default), parameter :: pi = 3.14159265358979323846d0
  integer, parameter :: log = 10

  ! methods
  public :: read_config
  public :: modify_config
  public :: print_config

contains

subroutine read_config
  print*, "config: reading ..."
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
  read(5,*) use_nwa ! 0 = actual top widths, 1 = tops in NWA
  read(5,*) sqrts
  read(5,*) seed
  read(5,*) itmx
  read(5,*) ncall
  read(5,*) nevents
  read(5,*) use_rambo
  read(5,*) map_phase_space
  read(5,*) symmetrise
  read(5,*) verbose
  read(5,*) ecm_low
  read(5,*) ecm_up
end subroutine read_config

subroutine modify_config

  integer :: i

  print*, "log:     ", trim(log_file)
  open(unit = log, file = log_file, status = "replace", action = "write")

  print*, "config: interpreting ..."

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
  print*, "config: printing to logfile ..."
  if (lhef_out) write(log,*) "events written to:", lhe_file
  if (ntuple_out) write(log,*) "events written to:", ntuple_file
  if (initial_state == 0) then
    if (final_state == -1) then
      write(log,*) "process: p p -> l+ l-"
    else if (final_state == 0) then
      write(log,*) 'process: p p -> t t~'
    else if (final_state == 1) then
      write(log,*) 'process: p p -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~'
    end if
  else if (initial_state == 1) then
    if (final_state == -1) then
      write(log,*) 'process: p p~ -> l+ l-'
    else if (final_state == 0) then
      write(log,*) 'process: p p~ -> t t~'
    else if (final_state == 1) then
      write(log,*) 'process: p p~ -> t t~ -> b b~ W+ W- -> b b~ l+ l- vl vl~'
    end if
  end if
  write(log,*) "preliminary vamp calls: ", ncall / 10
  write(log,*) "preliminary vamp iterations: ", itmx + 1
  write(log,*) "vamp calls: ", ncall
  write(log,*) "vamp iterations: ", itmx - 1
  write(log,*) "symmetrise parton momentum fraction: ", symmetrise
  write(log,*) "map phase space: ", map_phase_space
  write(log,*) 'random number seed: ', seed
  write(log,*) 'collider energy: ', sqrts
  write(log,*) "E_cm low: ", ecm_low
  write(log,*) "E_cm up:  ", ecm_up
  write(log,*) 'NWA: ', use_nwa
  write(log,*) 'RAMBO: ', use_rambo
  write(log,*) 'model: ', model_name
  write(log,*) "include gg: ", include_gg
  write(log,*) "include qq: ", include_qq
  write(log,*) "include uu: ", include_uu
  write(log,*) "include dd: ", include_dd
  write(log,*) "include a: ", include_a
  write(log,*) "include z: ", include_z
  write(log,*) "include x: ", include_x
  if (interference == 0) write(log,*) "interference: (gamma) + (Z) + (Z')"
  if (interference == 1) write(log,*) "interference: (gamma + Z) + (Z')"
  if (interference == 2) write(log,*) "interference: (gamma + Z + Z')"
  if (interference == 3) write(log,*) "interference: (gamma + Z + Z') - (gamma) - (Z)"
  if (interference == 4) write(log,*) "interference: (gamma + Z + Z') - (gamma) - (Z) - (Z')"
  if (ipdf ==  1) write(log,*) 'PDFs: CTEQ6m'
  if (ipdf ==  2) write(log,*) 'PDFs: CTEQ6d'
  if (ipdf ==  3) write(log,*) 'PDFs: CTEQ6l'
  if (ipdf ==  4) write(log,*) 'PDFs: CTEQ6l1'
  if (ipdf ==  5) write(log,*) 'PDFs: MRS99 (cor01)'
  if (ipdf ==  6) write(log,*) 'PDFs: MRS99 (cor02)'
  if (ipdf ==  7) write(log,*) 'PDFs: MRS99 (cor03)'
  if (ipdf ==  8) write(log,*) 'PDFs: MRS99 (cor04)'
  if (ipdf ==  9) write(log,*) 'PDFs: MRS99 (cor05)'
  if (ipdf == 10) write(log,*) 'PDFs: CT14LN'
  if (ipdf == 11) write(log,*) 'PDFs: CT14LL'
end subroutine print_config

end module configuration
