module configuration
  
  implicit double precision(a-h,o-z)
  
  common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  common/rndm/iseed
  common/reslocal/resl(20),standdevl(20)

  real    :: collider_energy
  integer :: initial_state ! 0 = pp, 1 = ppbar
  integer :: ifinal_state  ! 0 = no decay, 1 = dilepton, 2 = semi-had, 3 = full-had)
  integer :: nfinal
  integer :: o_structure
  character(50) :: model
  integer :: include_qcd
  integer :: include_ew
  integer :: include_bsm
  integer :: interference  
  integer :: o_nwa ! nwa flag (0:actual top widths,1: tops in nwa)
  !   branching ratio flag
  integer :: o_br
  !   transverse mass variables flag
  integer :: o_trans
  !   asymmetry observable flag
  integer :: o_asyms
  !   cut on top rapidity
  real :: ytmax
  !   cut on top pair boost
  real :: yttmin
  !   symmatrise of x1 and x2
  integer :: o_symx1x2
  !   symmatrise of x1 and x2
  integer :: o_symcost
  !   standard distributions flag
  integer :: o_distros
  !   2d-distributions flag
  integer :: o_dist2d
  !   set |m|^2=1
  integer :: o_m_eq_1

  ! distributions in asymmetries
  integer :: nasym
  parameter (nasym=9)
  integer :: nspat
  parameter (nspat=6) ! nasym-3
  integer :: o_asym(nasym)

  integer :: ixmax,jxmax
  
  public :: read_config
  public :: modify_config

  contains

    subroutine read_config

      ! read config file
      
      read(5,*) initial_state ! 0 = pp, 1 = ppbar
      read(5,*) collider_energy
      read(5,*) o_structure
      read(5,*) model
      read(5,*) include_qcd
      read(5,*) include_ew
      read(5,*) include_bsm
      read(5,*) interference

      read(5,*) ifinal_state ! 1=no decay,1=dilepton,2=semi-had,4=full-had)
      !   nwa flag ()
      read(5,*) o_nwa ! 0:actual top widths,1: tops in nwa
      !   branching ratio flag
      read(5,*) o_br
      !   transverse mass variables flag
      read(5,*) o_trans
      !   asymmetry observable flag
      read(5,*) o_asyms
      !   cut on top rapidity
      read(5,*) ytmax
      !   cut on top pair boost
      read(5,*) yttmin
      !   random number seed
      read(5,*) iseed
      !   maximum number of vegas iterations
      read(5,*) itmx
      !   number of vegas calls per iteration
      read(5,*) ncall
      !   desired accuracy (if negative, run maximum iterations.)
      read(5,*) acc
      !   symmatrise of x1 and x2
      read(5,*) o_symx1x2
      !   symmatrise of x1 and x2
      read(5,*) o_symcost
      !   standard distributions flag
      read(5,*) o_distros
      !   2d-distributions flag
      read(5,*) o_dist2d
      !   set |m|^2=1
      read(5,*) o_m_eq_1
    end subroutine read_config

    subroutine modify_config
      ! interpret config
      ! number of external lines
      if(ifinal_state == 0)then
        nfinal=4
      else
        nfinal=8
      end if
      ! nwa only for six-body final state
      if(ifinal_state == 0) o_nwa=0
      ! itmx no more than 20.
      if(itmx > 20)then
        write(*,*)'itmx does not have to exceed 20!'
        stop
      end if
      ! for every point in phase space with x1 and x2, include the point
      ! in phase space with x1<->x2
      if(o_symx1x2 == 1)then
        ixmax=2
      else
        ixmax=1
      end if
      ! in phase space with cost->-cost
      if(o_symcost == 1)then
        jxmax=2
      else
        jxmax=1
      end if
      ! do tops decay?
      if(ifinal_state == 0)then
        o_decay=0
      else
        o_decay=1
      end if
      ! in phase space with cost->-cost
      if(o_m_eq_1 == 1)then
        include_qcd=0
        include_ew=0
        include_bsm=0
      else
        jxmax=1
      end if
    end subroutine modify_config

end module configuration
