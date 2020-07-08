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
 
      ! initialize additional calculated parameters, maybe eft-modified
      module eftcouple
          use types
          implicit none

          public

          include 'ewcouple.f'
          include 'qcdcouple.f'
          include 'masses.f'

          real(dp), public, save, protected :: ecossin
          real(dp), public, save, protected :: eftgw, gb

          public :: eftcouple_init

          contains

          subroutine eftcouple_init()
              implicit none

              gb = sqrt(esq)/sqrt(1-xw)
              ecossin = sqrt(gw**2 + gb**2)

          end subroutine

      end module
