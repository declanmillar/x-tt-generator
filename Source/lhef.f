module lhef

  implicit none

  integer, parameter :: lhe = 20

  ! subroutines
  public :: lhe_init
  public :: lhe_beam
  public :: lhe_particle
  public :: lhe_close

contains


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
end subroutine lhe_init


subroutine lhe_beam(idbmup, ebmup, pdfgup, pdfsup)
  integer idbmup(2) ! ID of beam particle 1 and 2 according to the PDG
  real ebmup(2) ! energy in GeV of beam particles 1 and 2
  integer pdfgup(2) ! author group for beam 1 and 2 according to Cernlib PDFlib
  integer pdfsup(2) ! PDF set ID for beam 1 and 2 according to Cernlib PDFlib

  write(lhe,*) "<init>"
  write(lhe,*) idbmup, ebmup, pdfgup, pdfsup
  write(lhe,*) "</init>"
end subroutine lhe_beam


subroutine lhe_event()

end subroutine lhe_event


subroutine lhe_particle()

end subroutine lhe_particle


subroutine lhe_close
  write(lhe,*) "</LesHouchesEvents>"
end subroutine lhe_close

end module lhef
