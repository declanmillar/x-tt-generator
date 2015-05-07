module distributions

  use configuration
  use mathematics, only: pi
  use scattering, only: sigma
  use integration, only: it, cnorm
  use class_histogram
  use class_histogram2d

  implicit none

  public

! 1d distributions
  type(histogram) :: h_mtt
  type(histogram) :: h_mtt_reco
  type(histogram) :: h_mtb
  type(histogram) :: h_mt_reco
  type(histogram) :: h_beta
  type(histogram) :: h_cost
  type(histogram) :: h_et
  type(histogram) :: h_delta_y
  type(histogram) :: h_fl
  type(histogram) :: h_cosfl
  type(histogram) :: h_cost5
  type(histogram) :: h_cost7
  type(histogram) :: h_ct7ct5

  ! distribution in sigp
  real :: sigpmax,sigpmin,sigpw
  real :: xsigp(1000),fxsigp(n_asymmetries,1000,20) ,fxsigptot(n_asymmetries,1000)
  integer :: o_sigp
  integer :: ndiv_sigp

  ! distribution in sigm
  real :: sigmmax,sigmmin,sigmw
  real :: xsigm(1000),fxsigm(n_asymmetries,1000,20) ,fxsigmtot(n_asymmetries,1000)
  integer :: o_sigm
  integer :: ndiv_sigm

  public :: create_distributions
  public :: initialise_distributions
  public :: finalise_distributions

  integer, private :: i, j, k, iasy, jasy
  real, private :: sfxsigptot(n_asymmetries), sfxsigmtot(n_asymmetries), asym_int(n_asymmetries)

contains

subroutine create_distributions

  implicit none 

  print *, "Creating distributions..."


  if (o_mtt == 1) h_mtt = histogram("mtt", "d#sigma-/dM_{tt}--[pb]", "M_{tt}", 10, 0.d0, 13000.d0)
  if (o_mtb == 1) h_mtt = histogram("Mtt", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 140, 0.d0, 13000.d0)
  if (o_mt_reco == 1) h_mtt_reco = histogram("Mtt_reco", "d#sigma-/dM_{tt}--[pb/GeV]", 'M_{tt}--[GeV]', 140, 0.d0, 13000.d0)
  if (o_beta == 1) h_beta = histogram("beta", "d#sigma-/d#beta_{tt}--[pb/GeV]", 'beta_{t}', 100, 0.d0, 1000.d0)
  if (o_mtb == 1) h_mtb = histogram("mtb", "d#sigma-/d#m_{#bar{t}}--[pb]", "#m_{#bar{t}}", 100, 0.d0, 1000.d0)
  if (o_mtt_reco == 1) h_mt_reco = histogram("mt_reco", "d#sigma-/dm^{reco}_{t}--[pb]", "m^{reco}_{t}", 100, 0.d0, 1000.d0)
  if (o_beta == 1) h_beta = histogram("beta", "d#sigma-/d#beta--[pb]", "#beta", 100, 0.d0, 50.d0)
  if (o_cost == 1) h_cost = histogram("cost", "d#sigma-/dcos#theta_{t}--[pb]", "cos#theta_{t}", 100, -1.d0, 1.d0)
  if (o_et == 1) h_et = histogram("Et","d#sigma-/dE_{t}--[pb]", "E_{t}", 100, 0.d0, 1000.d0)
  if (o_delta_y == 1) h_delta_y = histogram("delta_y", "d#sigma-/d#Delta-y--[pb]", "#Delta-y", 100, -4.d0, 4.d0)
  if (o_fl == 1) h_fl = histogram("phil", "d#sigma-/d#phi_l--[pb]", "#phi_{l}", 100, 0.d0, 1000.d0)
  if (o_cosfl == 1) h_cosfl = histogram("cosphil", "d#sigma-/dcos#phi_{l}--[pb]", "cos#phi_{l}", 100, -1.d0, 1.d0)
  if (o_cost5 == 1) h_cost5 = histogram("cost5", "d#sigma-/dcos#theta_{l^{+}}--[pb]", "cos#theta_{l^{+}}", 100, -1.d0, 1.d0)
  if (o_cost7 == 1) h_cost7 = histogram("cost7", "d#sigma-/dcos#theta_{l^{-}}--[pb]", "cos#theta_{l^{-}}", 100, -1.d0, 1.d0)
  if (o_ct7ct5 == 1) h_ct7ct5 = histogram("ct7ct5","d#sigma-/dcos#theta_{l^{+}}cos#theta_{l^{-}}--[pb]", &
     "cos#theta_{l^+}cos#theta_{l^-}", 100,-1.d0,1.d0)

  ! sigp
  o_sigp = 1
  sigpmax = 13000
  sigpmin = 0
  ndiv_sigp = 100/5
  ! sigm
  o_sigm = 1
  sigmmax = 13000
  sigmmin = 0
  ndiv_sigm = 100/5

  do i = 1, n_asymmetries
    o_asym = 1
  end do

  if (final_state > 0) then
    o_asym(1) = 0
    o_asym(2) = 0
    o_asym(3) = 0

  else if (final_state == 0) then
    o_asym(6) = 0
    o_asym(10) = 0
    o_asym(11) = 0
    o_asym(12) = 0
  end if


  print *, "done."
end subroutine create_distributions

subroutine initialise_distributions

  implicit none

  print*, "Initialising distributions..."

  if (o_mtt == 1) call h_mtt%initialise()
  if (o_mtb == 1) call h_mtb%initialise()
  if (o_mt_reco == 1) call h_mt_reco%initialise()
  if (o_mtt_reco == 1) call h_mtt_reco%initialise()
  if (o_beta == 1) call h_beta%initialise()
  if (o_cost == 1) call h_cost%initialise()
  if (o_et == 1) call h_et%initialise()
  if (o_delta_y == 1) call h_delta_y%initialise()
  if (o_fl == 1) call h_fl%initialise()
  if (o_cosfl == 1) call h_cosfl%initialise()
  if (o_cost5 == 1) call h_cost5%initialise()
  if (o_cost7 == 1) call h_cost7%initialise()
  if (o_ct7ct5 == 1) call h_ct7ct5%initialise()

  if(o_sigp == 1)then
    sigpw=(sigpmax-sigpmin)/ndiv_sigp
    do i=1,ndiv_sigp
      xsigp(i)=sigpmin+sigpw*(i-1)+sigpw/2.d0
    end do
  end if

  if(o_sigm == 1)then
    sigmw=(sigmmax-sigmmin)/ndiv_sigm
    do i=1,ndiv_sigm
      xsigm(i)=sigmmin+sigmw*(i-1)+sigmw/2.d0
    end do
  end if

  print *, "done."

end subroutine initialise_distributions

subroutine finalise_distributions

  implicit none

  integer :: ndiv_sig

  print*, "Printing histograms..."
  open(unit = 10, file = 'Output/'//output_file, status = "replace", action = "write")

  write(10,*)'HISTOGRAMS'
    
  if (o_mtt == 1) call h_mtt % finalise()
  if (o_mtt_reco == 1) call h_mtt_reco % finalise()
  if (o_mtb == 1) call h_mtb % finalise()
  if (o_mt_reco == 1) call h_mt_reco % finalise()
  if (o_beta == 1) call h_beta % finalise()
  if (o_cost == 1) call h_cost % finalise()
  if (o_et == 1) call h_et % finalise()

  if (o_delta_y == 1) call h_delta_y % finalise()
  if (o_fl == 1) call h_fl%finalise()
  if (o_cosfl == 1) call h_cosfl%finalise()
  if (o_cost5 == 1) call h_cost5%finalise()
  if (o_cost7 == 1) call h_cost7%finalise()
  if (o_ct7ct5 == 1) call h_ct7ct5%finalise()

  ! plot distributions in all asymmetries
  if((o_sigp == 1) .and. (o_sigm == 1))then
    do jasy=1,n_asymmetries
      if(o_asym(jasy) == 0)then
        continue
      else
        ! snorm(jasy)=0.d0
        sfxsigptot(jasy)=0d0
        do j=1,ndiv_sigp
          fxsigptot(jasy,j)=0.d0
          do i=1,it
            fxsigp(jasy,j,i)=fxsigp(jasy,j,i)*sigma/cnorm(i)/sigpw
            fxsigptot(jasy,j)=fxsigptot(jasy,j)+fxsigp(jasy,j,i)
          end do
          sfxsigptot(jasy)=sfxsigptot(jasy)+fxsigptot(jasy,j)*sigpw
        end do
        sfxsigmtot(jasy)=0d0
        do j=1,ndiv_sigm
          do i=1,it
            fxsigm(jasy,j,i)=fxsigm(jasy,j,i)*sigma/cnorm(i)/sigmw
            fxsigmtot(jasy,j)=fxsigmtot(jasy,j)+fxsigm(jasy,j,i)
          end do
          sfxsigmtot(jasy)=sfxsigmtot(jasy)+fxsigmtot(jasy,j)*sigmw
        end do
        write(10,*)'ASYMMETRY'
        if(jasy == 1)then
          write(10,*)'ALL'
          write(10,*)'A_{ll}'
        else if(jasy == 2)then
          write(10,*)'AL'
          write(10,*)'A_{l}'
        else if(jasy == 3)then
          write(10,*)'APV'
          write(10,*)'A_{pv}'
        else if(jasy == 4)then
          write(10,*)'AFB'
          write(10,*)'A_{fb}'
        else if(jasy == 5)then
          write(10,*)'AFBstar'
          write(10,*)'A_{fb^{*}}'
        else if(jasy == 6)then
          write(10,*)'AFBstar_reco'
          write(10,*)'A_{fb^{*}}(reco)'
        else if(jasy == 7)then
          write(10,*)'AtRFB'
          write(10,*)'a^{t}_{rfb}'
        else if(jasy == 8)then
          write(10,*)'AttbRFB'
          write(10,*)"a^{b\bar{b}}_{rfb}"
        else if(jasy == 9)then
          write(10,*)'ARFB'
          write(10,*)"A_{rfb}"
         else if(jasy == 10)then
          write(10,*)'ARFB_reco'
          write(10,*)"A_{rfb}(reco)"
        else if(jasy == 11)then
          write(10,*)'A_l'
          write(10,*)'A_{l^+}'
        else if(jasy == 12)then
          write(10,*)'AlFB'
          write(10,*)'A^{l^{+}}_{FB}'
        end if
        write(10,*)'M_{tt}'
        ndiv_sig=(ndiv_sigm+ndiv_sigp)/2
        do i=1,ndiv_sig
          if(fxsigptot(jasy,i)+fxsigmtot(jasy,i) == 0.d0)then
            write(10,*)(xsigm(i)+xsigp(i))/2.d0,0.d0,0.d0,0.d0
            !           snorm(jasy)=snorm(jasy)+0.d0
          else
            write(10,*)(xsigm(i)+xsigp(i))/2.d0, &
            (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/ &
            (fxsigptot(jasy,i)+fxsigmtot(jasy,i)), &
            fxsigptot(jasy,i),fxsigmtot(jasy,i)
            !             snorm(jasy)=snorm(jasy)+
            !    &               (fxsigptot(jasy,i)-fxsigmtot(jasy,i))/
            !    &               (fxsigptot(jasy,i)+fxsigmtot(jasy,i))
            !    &               *fxmtttot(i)*mttw/sigma
          end if
        end do
        asym_int(jasy)=(sfxsigptot(jasy)-sfxsigmtot(jasy))/(sfxsigptot(jasy)+sfxsigmtot(jasy))
        write(10,*)'END'
        !         write(10,*)'(total asymmetry:',asym_int(jasy),')'
        !         write(10,*)'(integrated asymmetry:',snorm(jasy),' )'
      end if
    end do
  end if
  write(10,*) 'CLOSE'
  close(10)
  print *, "...complete."
end subroutine finalise_distributions

end module distributions
