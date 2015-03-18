program zprime

  ! calculates the cross section and generates distributions for
  ! pp -> tt,
  ! pp -> tt -> bw^+bbarw^- -> bbbare^+nue^-nubarc
  ! (future:) pp -> bw^+bbarw^- -> b bbar e^+ nu qqbar'
  ! uses adapted madgraph functions.
  ! uses cteq6 and mrs99 pdf subroutines.
  ! authors: declan millar, stefano moretti.

  use configuration
  use quantum_field_theory
  use kinematics
  use distributions
  use integration

  implicit none

  external differential_cross_section

  real :: o_width(5)
  real :: poltot(-1:1,-1:1),polchi(-1:1,-1:1)
  real :: spattot(n_fb_asymmetries,-1:1),spatchi(n_fb_asymmetries,-1:1)

  real :: avgi, chi2a, cross, error, sd, stantot, alfas, qcdl4
  integer :: iab, ndimensions, iasy, icteq, iphel,ispat,jphel

  
  ! branching ratio for t->bev=bmv=btv (with qcd corrections?)
  real :: brtbln=0.10779733d0
  ! branching ratio for t->bev=bmv=btv=1/9 (tree level)
  ! real brtbln/0.11111111d0/
  ! branching ratio for t->beq=bmq=btq=6/9 (tree level)
  ! real :: brtbeq=0.66666666d0

  integer :: i, j, k
  integer :: today(3), now(3)

  call read_config

  call modify_config

  call initialise_distributions

  s=collider_energy*collider_energy

  ! qcdl4 is qcd lambda4 (to match pdf fits).
  ! (pdfs are intrinsically linked to the value of lamda_qcd; alpha_qcd)
  if(structure_function == 1)qcdl4=0.326d0
  if(structure_function == 2)qcdl4=0.326d0
  if(structure_function == 3)qcdl4=0.326d0
  if(structure_function == 4)qcdl4=0.215d0
  if(structure_function == 5)qcdl4=0.300d0
  if(structure_function == 6)qcdl4=0.300d0
  if(structure_function == 7)qcdl4=0.300d0
  if(structure_function == 8)qcdl4=0.229d0
  if(structure_function == 9)qcdl4=0.383d0
  rlambdaqcd4=qcdl4
  ! initialise cteq grids.
  if(structure_function <= 4)then
    icteq=structure_function
    call setctq6(icteq)
  end if

  ! use appropriately evolved alphas.
  if(structure_function == 1)nloops=2
  if(structure_function == 2)nloops=2
  if(structure_function == 3)nloops=1
  if(structure_function == 4)nloops=1
  if(structure_function == 5)nloops=1
  if(structure_function == 6)nloops=1
  if(structure_function == 7)nloops=1
  if(structure_function == 8)nloops=1
  if(structure_function == 9)nloops=1

  ! initialise madgraph - masses and coupling constants of particles
  call initialise_standard_model
  call initialise_zprimes
  
  ! dimensions
  if(final_state == 0)then
    ndimensions=3
  else if(final_state > 0)then
    ndimensions=15
  end if
  ! if nprn<0 no print-out
  nprn=0

  if(final_state == 0)then
    m3=fmass(11)
    m4=m3
    m5=0.d0
    m6=0.d0
    m7=0.d0
    m8=0.d0
    ! integrates on:

    ! x(3) = (x1-tau)/(1-tau),
    ! x(2) = (ecm-m3-m4)/(ecm_max-m3-m4),
    ! x(1) = cos(theta3_cm)

    ! limits:
    do i=3,2,-1
      xl(i)=0.d0
      xu(i)=1.d0
    end do
    do i=1,1
      xl(i)=-1.d0
      xu(i)=1.d0
    end do

  else if(final_state >= 1)then
    m3=fmass(12)
    m4=m3
    m5=0.d0
    m6=0.d0
    m7=0.d0
    m8=0.d0
    ! integrates on:
         
    ! x(15) = (x1-tau)/(1-tau),
    ! x(14) = (ecm-m3-m4-m5-m6-m7-m8)
    !        /(ecm_max-m3-m4-m5-m6-m7-m8),
    ! x(13) = (xx356-xx356min)/(xx356max-xx356min),
    ! where xx356 = arctg((m356**2-m3**2)/m3/gamt),
    ! x(12) = (xx478-xx478min)/(xx478max-xx478min),
    ! where xx478 = arctg((m478**2-m3**2)/m3/gamt),
    ! x(11) = (xx56-xx56min)/(xx56max-xx56min),
    ! where xx56 = arctg((m56**2-rm_w**2)/rm_w/gamw),
    ! x(10) = (xx78-xx78min)/(xx78max-xx78min),
    ! where xx78 = arctg((m78**2-rm_w**2)/rm_w/gamw),
    ! x(9) = cos(theta_cm_356) = -cos(theta_cm_478)
    ! x(8) = cos(theta56_cm_356),
    ! x(7) = cos(theta78_cm_478),
    ! x(6) = cos(theta5_cm_56),
    ! x(5) = cos(theta7_cm_78),
    ! x(4) = fi56_cm_356,
    ! x(3) = fi78_cm_478,
    ! x(2) = fi5_cm_56,
    ! x(1) = fi8_cm_78;

    ! limits:
    do i=15,14,-1
      xl(i)=0.d0
      xu(i)=1.d0
    end do
    do i=13,10,-1
      xl(i)=0.d0
      xu(i)=1.d0
    end do
    do i=9,5,-1
      xl(i)=-1.d0
      xu(i)=1.d0
    end do
    do i=4,1,-1
      xl(i)=0.d0
      xu(i)=2.d0*pi
    end do
  end if

  call generate_bins

  ! output information before integration
  write(*,*)'====================================================='
  call idate(today)     ! today(1)=day, (2)=month, (3)=year
  call itime(now)       ! now(1)=hour, (2)=minute, (3)=second
  write(*,*)'date ',today(3),today(2),today(1)
  write(*,*)'time ',now(1),now(2),now(3)
  write(*,*)'-----------------------------------------------------'
  write(*,*)'process'
  if(initial_state == 0)then
    if(final_state == 0) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    ' #times br(t#rightarrow bl#nu)^{2}'
    if(final_state == 1) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    if(final_state == 2) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} q#bar{q} l #nu'
    if(final_state == 3) &
    write(*,*)'pp #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
  else if(initial_state == 1)then
    if(final_state == 0) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    ' #times br(t#rightarrow bl#nu)^{2}'
    if(final_state == 1) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} l^{+}l^{-} #nu#bar{#nu}'
    if(final_state == 2) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    '#rightarrow b#bar{b} q#bar{q} l #nu'
    if(final_state == 3) &
    write(*,*)'p#bar{p} #rightarrow t#bar{t}', &
    '#rightarrow b#bar{b} w^{+}w^{-}', &
    "#rightarrow b#bar{b} q#bar{q}q'#bar{q}'"
  end if
  write(*,*)'-----------------------------------------------------'
  write(*,*)'notes'
  write(*,*)'units: gev'
  write(*,*)'quarks: all massless except t, b.'
  if(structure_function == 1)write(*,*)'pdfs: cteq6m.'
  if(structure_function == 2)write(*,*)'pdfs: cteq6d.'
  if(structure_function == 3)write(*,*)'pdfs: cteq6l.'
  if(structure_function == 4)write(*,*)'pdfs: cteq6l1.'
  if(structure_function == 5)write(*,*)'pdfs: mrs99 (cor01).'
  if(structure_function == 6)write(*,*)'pdfs: mrs99 (cor02).'
  if(structure_function == 7)write(*,*)'pdfs: mrs99 (cor03).'
  if(structure_function == 8)write(*,*)'pdfs: mrs99 (cor04).'
  if(structure_function == 9)write(*,*)'pdfs: mrs99 (cor05).'
  if((final_state >= 1) .and. (use_nwa == 0))write(*,*)'tops: off-shell.'
  if((final_state >= 1) .and. (use_nwa == 1))write(*,*)'tops: nwa.'
  write(*,*)'bsm model_name: ',model_name
  if(include_qcd == 1)write(*,*)'qcd: on '
  if(include_qcd == 0)write(*,*)'qcd: off'
  if(include_ew == 1) write(*,*)'ew:  on '
  if(include_ew == 0) write(*,*)'ew:  off'
  if(include_bsm == 1)write(*,*)'bsm: on '
  if(include_bsm == 0)write(*,*)'bsm: off'
  if(interference == 0)write(*,*)'interference: none'
  if(interference == 1)write(*,*)'interference: sm'
  if(interference == 2)write(*,*)'interference: full'
  if(interference == 3)write(*,*)'interference: no square terms.'
  if(phase_space_only == 1)write(*,*)'phase space only'
  if(symmetrise_x1x2 == 1)write(*,*)'symmetrical integration over x1<->x2'
  if(symmetrise_costheta_t == 1)write(*,*)'symmetrical integration over cost'
  write(*,*)'seed: ',seed
  write(*,*)'-----------------------------------------------------'
  write(*,*)'parameters'
  write(*,*)'#sqrt{s}              ',collider_energy
  write(*,*)'at |y| <              ',abs(ytmax)
  write(*,*)'loops a_s evaluated at',nloops
  write(*,*)'a_{s}(m_{z})          ',alfas(rm_z,rlambdaqcd4,nloops)
  write(*,*)'#lambda_{qcd}(4)      ',qcdl4
  write(*,*)'m_{b}                 ',fmass(12)
  write(*,*)'#gamma_{b}            ',fwidth(12)
  write(*,*)'m_{t}                 ',fmass(11)
  write(*,*)'#gamma_{t}            ',fwidth(11)
  write(*,*)'m_{z}                 ',rm_z
  write(*,*)'#gamma_{z}            ',gamma_z
  write(*,*)'m_{w}                 ',rm_w
  write(*,*)'#gamma_{w}            ',gamma_w
  write(*,*)'m_{h}                 ',rm_h
  write(*,*)'#gamma_{h}            ',gamma_h
  write(*,*)'-----------------------------------------------------'
  write(*,*)'zprime parameters'
  do i=1,5
    if(rmzp(i) > 0)then
      write(*,*)'z#prime               ',i
      write(*,*)'m_{z#prime}           ',rmzp(i)
      write(*,*)'o_width:              ',o_width(i)
      write(*,*)'#gamma_{z#prime}      ',gamzp(i)
      write(*,*)'g_{p}                 ',gp(i)
      write(*,*)'g_{v}^{u}             ',gv_u(i)
      write(*,*)'g_{a}^{u}             ',ga_u(i)
      write(*,*)'g_{v}^{d}             ',gv_d(i)
      write(*,*)'g_{a}^{d}             ',ga_d(i)
      write(*,*)
    end if
  end do
  write(*,*)'-----------------------------------------------------'
  write(*,*)'cuts'
  write(*,*)'-----------------------------------------------------'


  ! reset counter
  npoints=0

  ! reset 
  if(final_state == 0)then
    do i=1,20
      resl(i)=0.d0
      standdevl(i)=0.d0
      cnorm(i)=0.d0
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          xsec_polar(i,iphel,jphel)=0.d0
          error_polar(i,iphel,jphel)=0.d0
        end do
      end do
      do ispat=1,n_fb_asymmetries
        do iasy=-1,+1,2
          xsec_fb(ispat,i,iasy)=0.d0
          error_fb(ispat,i,iasy)=0.d0
        end do
      end do
    end do
  end if

  ! integrate 
  write(*,*)'starting integration'
  it=0 
  call vegas(ndimensions, differential_cross_section, avgi, sd, chi2a)

  if(final_state == 0 .and. use_branching_ratio == 1)then
    ! multiply by branching ratios
    avgi=avgi*brtbln*brtbln
    sd=sd*brtbln*brtbln
  end if

  ! collect total cross-section
  cross=avgi
  sigma=avgi
  error=sd

  ! print integrated cross section
  write(*,*)''
  write(*,*)'integrated cross section'
  if(cross == 0d0)then
    write(*,*)'sigma = 0  ! check permitted gauge sectors.'
    stop
  else
    write(*,*)'sigma (pb)','error (same units)'
    write(*,*)cross,error
    write(*,*)'(using ',npoints,' points)'
  end if
  ! re-weight distributions for different iterations
  stantot=0.d0
  do i=1,it
    stantot=stantot+1.d0/standdevl(i)/standdevl(i)
  end do
  do i=1,it
    standdevl(i)=standdevl(i)*standdevl(i)*stantot
  end do
  do i=1,it
    cnorm(i)=resl(i)*standdevl(i)
  end do

  ! total asymmetries
  ! collect polarised cross sections.
  if(include_asymmetries == 1)then
    if(final_state == 0)then
      do iphel=-1,+1,2
        do jphel=-1,+1,2
          do i=1,it
            xsec_polar(i,iphel,jphel)=xsec_polar(i,iphel,jphel) &
            *avgi/cnorm(i)
            error_polar(i,iphel,jphel)=xsec_polar(i,iphel,jphel) &
            *sd/cnorm(i)
          end do
          poltot(iphel,jphel)=0.d0
          polchi(iphel,jphel)=0.d0
          do i=1,it
            poltot(iphel,jphel)=poltot(iphel,jphel) &
            +xsec_polar(i,iphel,jphel)
            polchi(iphel,jphel)=polchi(iphel,jphel) &
            +error_polar(i,iphel,jphel)
          end do
          polchi(iphel,jphel)=polchi(iphel,jphel) &
          /poltot(iphel,jphel)
          !        polchi(iphel,jphel)=
          ! & sqrt(abs(polchi(iphel,jphel)
          ! &         -poltot(iphel,jphel)**2*dfloat(ncall)))
          ! & /dfloat(ncall)
        end do
      end do
    end if

    ! collect unpolarised spatial asymmetry
    do ispat=1,n_fb_asymmetries
      if(o_asym(ispat+3) == 0)then
        continue
      else
        do iab=-1,+1,2
          do i=1,it
            xsec_fb(ispat,i,iab)=xsec_fb(ispat,i,iab) &
            *avgi/cnorm(i)
            error_fb(ispat,i,iab)=xsec_fb(ispat,i,iab) &
            *sd/cnorm(i)
          end do
          spattot(ispat,iab)=0.d0
          spatchi(ispat,iab)=0.d0
          do i=1,it     ! add up each iteration
            spattot(ispat,iab)=spattot(ispat,iab) &
            +xsec_fb(ispat,i,iab)
            spatchi(ispat,iab)=spatchi(ispat,iab) &
            +error_fb(ispat,i,iab)
          end do
          spatchi(ispat,iab)=spatchi(ispat,iab) &
          /spattot(ispat,iab)
          !         spatchi(iasy)=
          !    & sqrt(abs(spatchi(iasy)
          !    &         -spattot(iasy)**2*dfloat(ncall)))
          !    & /dfloat(ncall)
        end do
      end if
    end do

    ! define asymmetries
    if(final_state == 0)then
      ! all
      atot(1)= &
      +(poltot(+1,+1)-poltot(+1,-1) &
      -poltot(-1,+1)+poltot(-1,-1)) &
      /cross
      atoterr(1)= &
      +(polchi(+1,+1)+polchi(+1,-1) &
      +polchi(-1,+1)+polchi(-1,-1)) &
      /4.d0*atot(1)
      ! al
      atot(2)= &
      +(poltot(-1,-1)-poltot(+1,-1) &
      +poltot(-1,+1)-poltot(+1,+1)) &
      /cross
      atoterr(2)= &
      +(polchi(-1,-1)+polchi(+1,-1) &
      +polchi(-1,+1)+polchi(+1,+1)) &
      /4.d0*atot(2)
      ! apv
      atot(3)= &
      +(poltot(-1,-1)-poltot(+1,+1)) &
      /cross/2.d0
      atoterr(3)= &
      +(polchi(-1,-1)+polchi(+1,+1)) &
      /2.d0*atot(3)
    end if

    do iasy=4,n_asymmetries
      ispat=iasy-3
      if(o_asym(iasy) > 0)then
        atot(iasy)= &
        +(spattot(ispat,+1)-spattot(ispat,-1)) &
        /cross
        atoterr(iasy)= &
        +sd/avgi*atot(iasy)
      end if
    end do

    ! print asymmetries
    write(*,*)'total asymmetries'
    if(final_state == 0)then
      write(*,*)'ALL:                    uncertainty (same units):'
      write(*,*)atot(1),atoterr(1)
      write(*,*)'AL:                     uncertainty (same units):'
      write(*,*)atot(2),atoterr(2)
      write(*,*)'APV:                    uncertainty (same units):'
      write(*,*)atot(3),atoterr(3)
      write(*,*)'AFB:                    uncertainty (same units):'
      write(*,*)atot(4),atoterr(4)
      write(*,*)'AFB*:                   uncertainty (same units):'
      write(*,*)atot(5),atoterr(5)
      write(*,*)'AtRFB:                  uncertainty (same units):'
      write(*,*)atot(6),atoterr(6)
      write(*,*)"AttbRFB:                uncertainty (same units):"
      write(*,*)atot(7),atoterr(7)
      write(*,*)"ARFB/:                  uncertainty (same units):"
      write(*,*)atot(8),atoterr(8)
    else if(final_state > 0)then
      write(*,*)'A_l:                    uncertainty (same units):'
      write(*,*)atot(9),atoterr(9)
    end if
  end if

  call print_distributions
  call check_distributions

  stop
end program zprime
