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

      print *, "...done."
      
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
      
      print *, "...done."
    end subroutine modify_config

    subroutine setup_channels

      integer i, iasy

      print*, "Switching off irrelevent physics for chosen channel..."

      do iasy = 1, n_asymmetries
        o_asym(iasy) = include_asymmetries
      end do

      if (initial_state == 0) then
        ! disable non-useful variables for pp
        o_asym(4) = 0
        o_asym(7) = 0
        o_asym(8) = 0 
      end if

      if (initial_state == 1) then
        ! disable non-useful variables for ppbar
        o_asym(5) = 0
        o_asym(6) = 0
      end if

      if (final_state == 0) then
        ! disable 2to6 variables 
        o_ptb = 0
        o_ptbb = 0
        o_ptlp = 0
        o_ptlm = 0
        o_ptnu = 0
        o_ptnub = 0
        o_etab = 0
        o_etabb = 0
        o_etalp = 0
        o_etalm = 0
        o_etanu = 0
        o_etanub = 0
        o_phib = 0
        o_phibb = 0
        o_philp = 0
        o_philm = 0
        o_phinu = 0
        o_phinub = 0
        o_ycolb = 0
        o_ycolbb = 0
        o_ycollp = 0
        o_ycollm = 0
        o_ycolnu = 0
        o_ycolnub = 0
        o_etmiss = 0
        o_fl = 0
        o_dphi = 0
        o_cosfl = 0
        o_cost7 = 0
        o_cost5 = 0
        o_ct7ct5 = 0
        o_asym(6) = 0
        o_asym(10) = 0
        o_asym(11) = 0
        o_asym(12) = 0
        o_mtt_reco = 0
        o_mt_reco = 0
        o_mtb = 0
        o_ht = 0
        o_mttvis = 0
        o_mt1 = 0
        o_mt2 = 0
        o_mt3 = 0
        o_mct1 = 0
        o_mct2 = 0
        o_mct3 = 0
        o_mlt = 0
        o_mlct = 0
      end if

      if (final_state > 0) then
        ! disable non 2to6 variables
        o_asym(1) = 0
        o_asym(2) = 0
        o_asym(3) = 0
      end if

      if (final_state == 1) then 
        ! disable non-useful variables in dileptonic
        o_mll = 1
        o_mtt_reco = 0 
        o_mt_reco = 0
        o_mtb = 0
        do i = 4, 10
          o_asym(i) = 0
        end do
      end if

      if (final_state == 2) then
        ! disable non-useful variables in semi-leptonic
        o_cost7 = 0
        o_ct7ct5 = 0
        o_dphi = 0
        o_etmiss = 0
        o_ht = 0
        o_mttvis = 0
        o_mt1 = 0
        o_mt2 = 0
        o_mt3 = 0
        o_mct1 = 0
        o_mct2 = 0
        o_mct3 = 0
        o_mlt = 0
        o_mlct = 0
        o_asym(5) = 0
        o_asym(9) = 0
      end if

      if (final_state == 3) then
        ! disable non-useful variables in fully hadronic
        o_mtt_reco = 0
        o_cost5 = 0
        o_cost7 = 0
        o_ct7ct5 = 0
        o_dphi = 0
        o_etmiss = 0
        o_mt_reco = 0
        o_mtb = 0
        o_ht = 0
        o_mttvis = 0
        o_mt1 = 0
        o_mt2 = 0
        o_mt3 = 0
        o_mct1 = 0
        o_mct2 = 0
        o_mct3 = 0
        o_mlt = 0
        o_mlct = 0

        o_asym(6) = 0
        o_asym(10) = 0
        o_asym(11) = 0
        o_asym(12) = 0
      end if

      ! disable A_PV
      o_asym(3) = 0

      print *, "...done."
      
    end subroutine setup_channels

end module configuration
