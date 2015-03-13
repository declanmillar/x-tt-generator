module Configuration
  
  implicit double precision(a-h,o-z)
  
  common/bveg1/ncall,itmx,nprn,ndev,xl(100),xu(100),acc
  common/bveg2/it,ndo,si,swgt,schi,xi(50,100)
  common/rndm/iseed
  common/reslocal/resl(20),standdevl(20)

  !   Collider flag (initial_state=0: pp; initial_state=1: ppbar)
  integer :: initial_state
  !   Collider energy
  real :: collider_energy
  !   PDFs
  integer :: o_structure
  !   Name of model file
  character(50) :: model
  !   Permitted gauge sector options
  integer :: include_QCD
  integer :: include_EW
  integer :: include_BSM
  !   Interference options
  integer :: interference
  !   Final state option (1:no decay,1:dilepton,2:semi-had,4:full-had)
  integer :: ifinal_state
  !   NWA flag (0:Actual top widths,1: tops in NWA)
  integer :: nfinal
  integer :: o_NWA
  !   Branching ratio flag
  integer :: o_BR
  !   Transverse mass variables flag
  integer :: o_trans
  !   Asymmetry observable flag
  integer :: o_asyms
  !   Cut on top rapidity
  real :: ytmax
  !   Cut on top pair boost
  real :: yttmin
  !   Symmatrise of x1 and x2
  integer :: o_symx1x2
  !   Symmatrise of x1 and x2
  integer :: o_symcost
  !   Standard distributions flag
  integer :: o_distros
  !   2d-distributions flag
  integer :: o_dist2d
  !   set |M|^2=1
  integer :: o_M_eq_1

  integer :: ixmax,jxmax
  
  public :: read_config
  public :: modify_config

  contains

    subroutine read_config

      ! Read config file
      
      read(5,*) initial_state ! 0 = pp, 1 = ppbar
      read(5,*) collider_energy
      read(5,*) o_structure
      read(5,*) model
      read(5,*) include_QCD
      read(5,*) include_EW
      read(5,*) include_BSM
      read(5,*) interference

      read(5,*) ifinal_state ! 1=no decay,1=dilepton,2=semi-had,4=full-had)
      !   NWA flag ()
      read(5,*) o_NWA ! 0:Actual top widths,1: tops in NWA
      !   Branching ratio flag
      read(5,*) o_BR
      !   Transverse mass variables flag
      read(5,*) o_trans
      !   Asymmetry observable flag
      read(5,*) o_asyms
      !   Cut on top rapidity
      read(5,*) ytmax
      !   Cut on top pair boost
      read(5,*) yttmin
      !   Random number seed
      read(5,*) iseed
      !   Maximum number of Vegas iterations
      read(5,*) itmx
      !   Number of Vegas calls per iteration
      read(5,*) ncall
      !   Desired Accuracy (If negative, run maximum iterations.)
      read(5,*) acc
      !   Symmatrise of x1 and x2
      read(5,*) o_symx1x2
      !   Symmatrise of x1 and x2
      read(5,*) o_symcost
      !   Standard distributions flag
      read(5,*) o_distros
      !   2d-distributions flag
      read(5,*) o_dist2d
      !   set |M|^2=1
      read(5,*) o_M_eq_1
    end subroutine read_config

    subroutine modify_config
      ! Interpret config
      ! Number of external lines
      if(ifinal_state == 0)then
        nfinal=4
      else
        nfinal=8
      end if
      ! NWA only for six-body final state
      if(ifinal_state == 0) o_NWA=0
      ! itmx no more than 20.
      if(itmx > 20)then
        write(*,*)'itmx does not have to exceed 20!'
        stop
      end if
      ! For every point in phase space with x1 and x2, include the point
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
      ! Do tops decay?
      if(ifinal_state == 0)then
        o_decay=0
      else
        o_decay=1
      end if
      ! in phase space with cost->-cost
      if(o_M_eq_1 == 1)then
        include_QCD=0
        include_EW=0
        include_BSM=0
      else
        jxmax=1
      end if
    end subroutine modify_config

end module Configuration
