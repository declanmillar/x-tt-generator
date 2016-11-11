! vamp_kinds.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module kinds
    implicit none
    integer, parameter, private :: single = &
        & selected_real_kind (precision(1.0), range(1.0))
    integer, parameter, private :: double = &
        & selected_real_kind (precision(1.0_single) + 1, range(1.0_single) + 1)
    integer, parameter, private :: extended = &
        & selected_real_kind (precision (1.0_double) + 1, range (1.0_double))
    integer, parameter, public :: default = double
    character(len=*), public, parameter :: KINDS_RCS_ID = &
        "$Id: kinds.nw 314 2010-04-17 20:32:33Z ohl $"
end module kinds