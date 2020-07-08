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
 
      module scaleset_m
          use types
          implicit none

          public :: set_scales

          contains

          subroutine set_scales(renscale_in,facscale_in)
              use constants, only: pi
              implicit none

              include 'couple.f'
              include 'qcdcouple.f'
              include 'nlooprun.f'
              include 'scale.f'
              include 'facscale.f'

              real(dp), intent(in) :: renscale_in, facscale_in

              real(dp) :: alphas

              scale = renscale_in
              facscale = facscale_in

              as = alphas(scale,amz,nlooprun)

              ason2pi = as/2/pi
              ason4pi = as/4/pi
              gsq = 4*pi*as
              musq = scale**2
          end subroutine
      end module
