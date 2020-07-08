!  Copyright (C) 2019, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
 
      module cxx11random
        implicit none

        interface
          subroutine cxx11_init_random(seeds) bind(C, name="cxx11_init_random")
            use ISO_C_BINDING, only: c_int, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: seeds
          end subroutine

          function cxx11_random_number() bind(C, name="cxx11_random_number")
            use ISO_C_BINDING, only: c_double
            implicit none
            real(c_double) :: cxx11_random_number
          end function
        end interface

      end module
