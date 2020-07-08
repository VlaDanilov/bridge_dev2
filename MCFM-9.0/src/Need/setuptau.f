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
 
      subroutine setuptau()
      implicit none
      include 'types.f'
      include 'taucut.f'
      include 'nproc.f'

      select case (nproc)
          case (1,6)
              taucut = 6.e-3_dp
          case (31,32)
              taucut = 6.e-3_dp
          case (91:110)
              taucut = 3.e-3_dp
          case (111,112,119)
              taucut = 4.e-3_dp
          case (285)
              taucut = 1.e-4_dp
          case (300,305)
              taucut = 3.e-4_dp
          case default
              write(6,*) 'Pre-determined value of taucut not available'
              write(6,*) 'for nproc = ',nproc
              stop
      end select

      end
      
