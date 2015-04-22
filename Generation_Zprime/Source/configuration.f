module configuration

  use integration, only: seed, itmx, ncall, acc

  implicit none

  ! parameters in config file
  real    :: collider_energy
  integer :: initial_state ! 0 = pp, 1 = ppbar
  integer :: final_state  ! 0 = no decay, 1 = dilepton, 2 = semi-had, 3 = full-had
  
  integer :: structure_function
  character(50) :: model_name
  character(100) :: ntuple_file
  character(50) :: output_file
  integer :: include_qcd
  integer :: include_ew
  integer :: include_bsm
  integer :: interference  
  integer :: use_nwa
  integer :: use_branching_ratio ! 0 = no, 1 = dilepton, 2 = semi-had, 3 = full-had
  integer :: include_transverse
  integer :: include_asymmetries
  real :: ytmax
  real :: yttmin
  integer :: symmetrise_x1x2
  integer :: symmetrise_costheta_t
  integer :: print_all_distributions
  integer :: print_2d_distributions
  integer :: include_errors
  integer :: phase_space_only

  ! derived parameters
  integer :: n_final
  integer :: tops_decay

  ! distributions in asymmetries
  integer n_asymmetries
  parameter (n_asymmetries = 12)
  integer n_fb_asymmetries
  parameter (n_fb_asymmetries = 9)
  integer :: o_asym(n_asymmetries)

  integer :: ixmax,jxmax
  
  public :: read_config
  public :: modify_config

  contains

    subroutine read_config

      ! read config file
      print *, "Reading config file..."

      read(5,*) ntuple_file

      read(5,*) output_file
      
      read(5,*) initial_state ! 0 = pp, 1 = ppbar

      read(5,*) final_state ! 1=no decay,1=dilepton,2=semi-had,4=full-had)
      
      read(5,*) model_name
      
      read(5,*) structure_function

      read(5,*) include_qcd

      read(5,*) include_ew

      read(5,*) include_bsm

      read(5,*) phase_space_only

      read(5,*) interference

      read(5,*) use_branching_ratio

      read(5,*) use_nwa ! 0:actual top widths,1: tops in nwa

      read(5,*) include_transverse

      read(5,*) include_asymmetries

      read(5,*) collider_energy

      read(5,*) ytmax

      read(5,*) yttmin

      read(5,*) seed

      read(5,*) itmx

      read(5,*) ncall

      read(5,*) acc

      read(5,*) symmetrise_x1x2

      read(5,*) symmetrise_costheta_t

      read(5,*) print_all_distributions

      read(5,*) print_2d_distributions

      read(5,*) include_errors

      print *, "done."
      
    end subroutine read_config

    subroutine modify_config

      print*, "Modifying config file..."

      if(final_state == 0)then
        n_final=4
      else
        n_final=8
      end if

      if(final_state == 0) use_nwa = 0
      if(final_state > 0) use_branching_ratio = 0

      if (itmx > 20 )then
        write(*,*) 'itmx does not have to exceed 20!'
        stop
      end if

      if(symmetrise_x1x2 == 1)then
        ixmax=2
      else
        ixmax=1
      end if

      if(symmetrise_costheta_t == 1)then
        jxmax=2
      else
        jxmax=1
      end if

      if(final_state == 0)then
        tops_decay=0
      else
        tops_decay=1
      end if

      if(phase_space_only == 1)then
        include_qcd=0
        include_ew=0
        include_bsm=0
      end if
      
      print *, "done."
    end subroutine modify_config

end module configuration
