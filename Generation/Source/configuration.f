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
  integer :: additional_kinematics
  integer :: symmetrise_x1x2
  integer :: symmetrise_costheta_t
  integer :: print_distributions
  integer :: include_errors
  integer :: phase_space_only
  integer :: use_rambo
  integer :: map_phase_space
  integer :: verbose
  integer :: nfb
  real :: lambdaqcd4
  integer :: nloops

  ! derived parameters
  integer :: n_final
  integer :: tops_decay

  ! distributions in asymmetries
  integer, parameter :: n_asymmetries = 12
  integer :: o_asym(n_asymmetries)

  integer :: ixmax, jxmax

  ! variable switches
  integer :: o_ptb = 1
  integer :: o_ptbb = 1
  integer :: o_ptlp = 1
  integer :: o_ptlm = 1
  integer :: o_ptnu = 1
  integer :: o_ptnub = 1
  integer :: o_ptt = 1
  integer :: o_pttb = 1
  integer :: o_etab = 1
  integer :: o_etabb = 1
  integer :: o_etalp = 1
  integer :: o_etalm = 1
  integer :: o_etanu = 1
  integer :: o_etanub = 1
  integer :: o_etat = 1
  integer :: o_etatb = 1
  integer :: o_phib = 1
  integer :: o_phibb = 1
  integer :: o_philp = 1
  integer :: o_philm = 1
  integer :: o_phinu = 1
  integer :: o_phinub = 1
  integer :: o_phit = 1
  integer :: o_phitb = 1
  integer :: o_ycolb = 1
  integer :: o_ycolbb = 1
  integer :: o_ycollp = 1
  integer :: o_ycollm = 1
  integer :: o_ycolnu = 1
  integer :: o_ycolnub = 1
  integer :: o_ycolt   = 1
  integer :: o_ycoltb = 1
  integer :: o_mtt = 1
  integer :: o_mtt_reco = 1
  integer :: o_mtb = 1
  integer :: o_mt_reco = 1
  integer :: o_etmiss = 1
  integer :: o_beta = 1
  integer :: o_cost = 1
  integer :: o_et = 1
  integer :: o_delta_y = 1
  integer :: o_fl = 1
  integer :: o_cosfl = 1
  integer :: o_dphi = 1
  integer :: o_cost5 = 1
  integer :: o_cost7 = 1
  integer :: o_ct7ct5 = 1
  integer :: o_mll = 0
  integer :: o_ht = 1
  integer :: o_mttvis = 1
  integer :: o_mt1 = 1
  integer :: o_mt2 = 1
  integer :: o_mt3 = 1
  integer :: o_mct1 = 1
  integer :: o_mct2 = 1
  integer :: o_mct3 = 1
  integer :: o_mlt = 1
  integer :: o_mlct = 1
  
  public :: read_config
  public :: modify_config

  contains

    subroutine read_config

      ! read config file
      print*, "Reading config file..."

      read(5,"(a)") ntuple_file

      read(5,*) output_file
      
      read(5,*) initial_state ! 0 = pp, 1 = ppbar

      read(5,*) final_state ! 1 = no decay, 1 = dilepton, 2 = semilepton, 4 = full hadron
      
      read(5,*) model_name
      
      read(5,*) structure_function

      read(5,*) include_qcd

      read(5,*) include_ew

      read(5,*) include_bsm

      read(5,*) phase_space_only

      read(5,*) interference

      read(5,*) use_nwa ! 0:actual top widths,1: tops in nwa

      read(5,*) additional_kinematics

      read(5,*) collider_energy

      read(5,*) seed

      read(5,*) itmx

      read(5,*) ncall

      read(5,*) acc

      read(5,*) use_rambo

      read(5,*) map_phase_space

      read(5,*) symmetrise_x1x2

      read(5,*) symmetrise_costheta_t

      read(5,*) print_distributions

      read(5,*) include_errors

      read(5,*) verbose

      print*, "...complete."
      
    end subroutine read_config

    subroutine modify_config

      print*, "Modifying config file..."

      if(final_state == 0) then
        nfb = 5
      else if (final_state > 0) then
        nfb = 9
      end if

      if(final_state == 0)then
        n_final = 4
      else
        n_final = 8
      end if

      if (use_rambo == 1) then
        map_phase_space = 0
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
      
      print*, "...complete."
    end subroutine modify_config

end module configuration
