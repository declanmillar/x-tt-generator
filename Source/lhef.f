module lhef

  implicit none

  integer, parameter :: lhe = 20

  ! subroutines
  public :: lhe_init
  public :: lhe_beam
  public :: lhe_process
  public :: lhe_add_event
  public :: lhe_add_particle
  public :: lhe_close

  private :: mass

contains

function mass(p)
  real :: mass, p(4)
  mass = sqrt(abs(p(4)*p(4) - p(1)*p(1) - p(2)*p(2) - p(3)*p(3)))
end function mass

subroutine lhe_init(lhe_file)
  character(*) lhe_file
  open(unit = lhe, file = lhe_file, status = "replace", action = "write")
  write(lhe,*) '<LesHouchesEvents version="1.0">'
  write(lhe,*) "<!--"
  write(lhe,*) "# File generated with zprime-top-generator"
  write(lhe,*) "-->"
  write(lhe,*) "<header>"
  write(lhe,*) "<!-- individually designed XML tags, in fancy XML style -->"
  write(lhe,*) "</header>"
  write(lhe,*) "<init>"
end subroutine lhe_init

subroutine lhe_beam(idbmup, ebmup, pdfgup, pdfsup)
  integer :: idbmup(2) ! ID of beam particle 1 and 2 according to the PDG
  real :: ebmup(2) ! energy in GeV of beam particles 1 and 2
  integer :: pdfgup(2) ! author group for beam 1 and 2 according to Cernlib PDFlib
  integer :: pdfsup(2) ! PDF set ID for beam 1 and 2 according to Cernlib PDFlib

  write(lhe,*) idbmup, ebmup, pdfgup, pdfsup
end subroutine lhe_beam

subroutine lhe_process(idwtup, nprup, xsecup, xerrup, xmaxup, lprup)
  integer :: idwtup ! master switch dictating how the event weights (XWGTUP) are interpreted
  integer :: nprup ! the number of different user subprocesses
  real :: xsecup(nprup) ! the cross section for process j in pb this entry is mandatory for idwtup=±2.
  real :: xerrup(nprup) ! the statistical error associated with the cross section of process j in pb
  real :: xmaxup(nprup) ! the maximum xwgtup for process j
  integer :: lprup(nprup) ! a listing of all user process IDs that can appear in IDPRUP of HEPEUP for this run

  write(lhe,*) idwtup, nprup, xsecup, xerrup, xmaxup, lprup
  write(lhe,*) "</init>"
end subroutine lhe_process

subroutine lhe_add_event(nup, idprup, xwgtup, scalup, aqedup, aqcdup)
  integer :: nup ! number of particle entries in this event
  integer :: idprup ! ID of the process for this event
  real :: xwgtup ! event weight
  real :: scalup ! scale of the event in GeV, as used for calculation of PDFs
  real :: aqedup ! the QED coupling αQED used for this event (e.g. 1/128)
  real :: aqcdup ! the QCD coupling αQCD used for this event

  write(lhe,*) "<event>"
  write(lhe,*) nup, idprup, xwgtup, scalup, aqedup, aqcdup
end subroutine lhe_add_event

subroutine lhe_add_particle(idup, istup, mothup1, mothup2, icolup1, icolup2, pup, vtimup, spinup)
  integer :: idup ! particle ID according to Particle Data Group convention
  integer :: istup ! status code
  integer :: mothup1, mothup2 ! index of first and last mother
  integer :: icolup1, icolup2 ! 1(2) is the integer tag for the (anti) color flow line passing through the color of the particle
  real :: pup(4) ! lab frame momentum
  real :: vtimup ! invariant lifetime cτ (distance from production to decay) in mm
  real :: spinup ! cosine of the angle between the spin-vector of particle I and the 3-momentum of the decaying particle, specified in the lab frame

  write(lhe,*) idup, istup, mothup1, mothup2, icolup1, icolup2, pup, mass(pup), vtimup, spinup
end subroutine lhe_add_particle

subroutine lhe_end_event
  write(lhe,*) "</event>"
end subroutine lhe_end_event


subroutine lhe_close
  write(lhe,*) "</LesHouchesEvents>"
end subroutine lhe_close

end module lhef
