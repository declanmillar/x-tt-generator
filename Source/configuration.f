module configuration

  use integration, only: seed, itmx, ncall, acc

  implicit none

  ! Config
  integer :: ntuple_out
  integer :: lhef_out
  real :: collider_energy
  integer :: initial_state
  integer :: final_state
  integer :: structure_function
  character(50) :: model_name
  character(100) :: ntuple_file
  character(100) :: log_file
  character(100) :: lhe_file
  integer :: include_signal
  integer :: include_background
  integer :: include_gg
  integer :: include_qq
  integer :: include_uu
  integer :: include_dd
  integer :: include_a
  integer :: include_z
  integer :: include_x
  integer :: interference
  integer :: use_nwa
  integer :: symmetrise
  integer :: phase_space_only
  integer :: use_rambo
  integer :: map_phase_space
  integer :: verbose
  real :: ecm_low, ecm_up
  integer :: cut

  ! Derived config
  integer :: n_final
  integer :: tops_decay
  integer :: ffinal
  integer :: nloops
  real :: lambdaqcd4
  integer :: z_mixing
  real :: s2mix

  ! constants
  real, parameter :: pi = 3.14159265358979323846d0
  integer, parameter :: log = 10

  ! Methods
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
    read(5,*) structure_function
    read(5,*) include_signal
    read(5,*) include_background
    read(5,*) include_gg
    read(5,*) include_qq
    read(5,*) include_uu
    read(5,*) include_dd
    read(5,*) include_a
    read(5,*) include_z
    read(5,*) include_x
    read(5,*) phase_space_only
    read(5,*) interference
    read(5,*) use_nwa ! 0:actual top widths,1: tops in nwa
    read(5,*) collider_energy
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

    i = len(log_file)
    do while(log_file(i:i) == '')
      i = i - 1
    end do
    log_file = log_file(1:i)
    print*, "Log:    ", log_file

    i = len(lhe_file)
    do while(lhe_file(i:i) == '')
      i = i - 1
    end do
    lhe_file = lhe_file(1:i)
    print*, "LHEF:   ", lhe_file

    i = len(ntuple_file)
    do while(ntuple_file(i:i) == '')
      i = i - 1
    end do
    ntuple_file = ntuple_file(1:i)
    print*, "Ntuple: ", ntuple_file


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
  if (structure_function == 1) write(log,*) 'PDFs:CTEQ6m'
  if (structure_function == 2) write(log,*) 'PDFs:CTEQ6d'
  if (structure_function == 3) write(log,*) 'PDFs:CTEQ6l'
  if (structure_function == 4) write(log,*) 'PDFs:CTEQ6l1'
  if (structure_function == 5) write(log,*) 'PDFs:MRS99 (cor01)'
  if (structure_function == 6) write(log,*) 'PDFs:MRS99 (cor02)'
  if (structure_function == 7) write(log,*) 'PDFs:MRS99 (cor03)'
  if (structure_function == 8) write(log,*) 'PDFs:MRS99 (cor04)'
  if (structure_function == 9) write(log,*) 'PDFs:MRS99 (cor05)'
  if (structure_function == 10) write(log,*) 'PDFs:CT14LN'
  if (structure_function == 11) write(log,*) 'PDFs:CT14LL'
  if ((final_state >= 1) .and. (use_nwa == 1)) write(log,*) 'NWA:ON'
  if ((final_state >= 1) .and. (use_nwa == 0)) write(log,*) 'NWA:OFF'
  write(log,*) 'Model:', model_name
  if (include_a == 1) write(log,*) 'photon:ON '
  if (include_a == 0) write(log,*) 'photon:OFF'
  if (include_z == 1) write(log,*) 'Z boson:ON '
  if (include_z == 0) write(log,*) 'Z boson:OFF'
  if (include_x == 1) write(log,*) 'BSM:ON '
  if (include_x == 0) write(log,*) 'BSM:OFF'
  write(log,*) "include_gg: ", include_gg
  write(log,*) "include_qq: ", include_qq
  write(log,*) "include_uu: ", include_uu
  write(log,*) "include_dd: ", include_dd
  if (interference == 0) write(log,*) "Interference:none"
  if (interference == 1) write(log,*) "Interference:(gamma + Z) + (Z')"
  if (interference == 2) write(log,*) "Interference:(gamma + Z + Z')"
  if (interference == 3) write(log,*) "Interference:(gamma + Z + Z') - (gamma) - (Z)"
  if (interference == 4) write(log,*) "Interference:(gamma + Z + Z') - (gamma) - (Z) - (Z')"
  if (symmetrise == 0) write(log,*) 'Not symmetrising integration: x1<->x2!'
  if (use_rambo == 1) write(log,*) 'RAMBO:ON'
  if (map_phase_space == 0) write(log,*) "Phase space mapping:ON"
  write(log,*) 'Seed:', seed
  write(log,*) 'Collider energy:', collider_energy
  if (ecm_low > 0) write(log,*) "E_CM low         :", ecm_low
  if (ecm_up > 0) write(log,*) "E_CM up             ", ecm_up
  if (lhef_out == 1) write(log,*) "Events written to:", lhe_file
  if (ntuple_out == 1) write(log,*) "Events written to:", ntuple_file
end subroutine print_config

end module configuration
