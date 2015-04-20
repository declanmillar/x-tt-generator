module class_Histogram2d
  use configuration, only: include_errors
  use kinematics, only: sigma
  use integration, only: it, cnorm
  implicit none
  private

  real :: diff_max = 1e-12

  type, public :: Histogram2d
    private
    character (len = 50), public :: name
    character (len = 50), public :: ztitle
    character (len = 50), public :: xtitle
    real, public :: xlow
    real, public :: xup
    integer, public :: nxbin
    character (len = 50), public :: ytitle
    real, public :: ylow
    real, public :: yup
    integer, public :: nybin

    real :: xbinw = 0, ybinw = 0
    real :: x(500) = 0, y(500) = 0
    real :: fxy(500,500,20) = 0, fxytot(500,500) = 0
    real :: sfxytot=0
    real :: sumw2(500,500,20) = 0, sumw2tot(500,500) = 0

    contains
      procedure :: bin_width => histogram2d_bin_width 
      procedure :: midpoints => histogram2d_midpoints
      procedure :: initialise => histogram2d_initialise
      procedure :: finalise => histogram2d_finalise
      procedure :: fill => histogram2d_fill
      procedure :: collate => histogram2d_collate
      procedure :: check => histogram2d_check
      procedure :: print => histogram2d_print
  end type Histogram2d

contains

  subroutine histogram2d_initialise(this)
    class(histogram2d), intent(inout) :: this
    call this%bin_width
    call this%midpoints
  end subroutine histogram2d_initialise

  subroutine histogram2d_finalise(this)
    class(histogram2d), intent(inout) :: this
    call this%collate
    call this%check
    call this%print
  end subroutine histogram2d_finalise

  subroutine histogram2d_bin_width(this)
    class(histogram2d), intent(inout) :: this
    this%xbinw = (this%xup-this%xlow)/this%nxbin
    this%ybinw = (this%yup-this%ylow)/this%nybin
  end subroutine histogram2d_bin_width

  subroutine histogram2d_midpoints(this)
    class(histogram2d), intent(inout) :: this
    integer :: i
    do i = 1, this%nxbin
      this%x(i) = this%xlow + this%xbinw * (i - 1) + this%xbinw / 2.d0
    end do
    do i = 1, this%nybin
      this%y(i) = this%ylow + this%ybinw * (i - 1) + this%ybinw / 2.d0
    end do
  end subroutine histogram2d_midpoints

  subroutine histogram2d_fill(this, xvalue, yvalue, weight)
    class(histogram2d), intent(inout) :: this
    real :: xvalue, yvalue, weight
    integer :: xbin, ybin

    xbin = int((xvalue - this%xlow)/this%xbinw) + 1
    ybin = int((yvalue - this%ylow)/this%ybinw) + 1

    if ((xbin > this%nxbin) .or. (ybin > this%nybin)) then
      continue
    else if ((xbin < 1) .or. (ybin < 1)) then
      continue
    else
      this%fxy(xbin, ybin, it) = this%fxy(xbin, ybin, it) + weight
      if (include_errors == 1) this%sumw2(xbin, ybin, it) = this%sumw2(xbin, ybin, it) + weight*weight
    end if

  end subroutine histogram2d_fill

  subroutine histogram2d_collate(this)

    ! average over iterations, accounting for weight 

    class(histogram2d), intent(inout) :: this
    integer :: i, j, k

    do i = 1, this%nxbin
      do j = 1, this%nybin
        do k = 1, it
          this%fxytot(i,j) = 0.d0
          this%fxy(i,j,k) = this%fxy(i,j,k)*sigma/cnorm(k)/this%xbinw/this%ybinw
          this%fxytot(i,j) = this%fxytot(i,j)+this%fxy(i,j,k)
          if (include_errors == 1) then 
            this%sumw2(i,j,k) = this%sumw2(i,j,k)*sigma/cnorm(k)/this%xbinw/this%ybinw*sigma/cnorm(k)/this%xbinw/this%ybinw
            this%sumw2tot(i,j) = this%sumw2tot(i,j)+this%sumw2(i,j,k)
          else
            this%sumw2tot(i,j) = 0
          end if
        end do
        this%sfxytot = this%sfxytot + this%fxytot(i,j)*this%xbinw*this%ybinw
      end do
    end do
  end subroutine histogram2d_collate

  subroutine histogram2d_check(this)
    class(histogram2d), intent(inout) :: this
      if(abs(sigma-this%sfxytot) > diff_max) then
        write(10,*) this%name, ' error: ', this%sfxytot
      end if
  end subroutine histogram2d_check

  subroutine histogram2d_print(this)
    class(histogram2d), intent(inout) :: this
    integer :: i, j
    write(10,*)'DISTRIBUTION'
    write(10,*) this%name
    write(10,*) this%ztitle
    write(10,*) this%xtitle
    write(10,*) this%xlow
    write(10,*) this%xup
    write(10,*) this%nxbin
    write(10,*) this%ytitle
    write(10,*) this%ylow
    write(10,*) this%yup
    write(10,*) this%nybin
    do i = 1, this%nxbin
        do j = 1, this%nybin
          write(10,*) this%x(i), this%y(j), this%fxytot(i,j)
        end do
      end do
    write(10,*)'END'
  end subroutine histogram2d_print

end module class_Histogram2d
