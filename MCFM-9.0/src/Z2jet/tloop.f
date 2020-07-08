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
 
      function tloop(s23,mtsq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      complex(dp):: tloop
      real(dp):: s23,mtsq
      complex(dp):: Bdiff

      tloop=-1._dp-2._dp*mtsq/s23*(
     & 3._dp*s23*qlI3(zip,zip,s23,mtsq,mtsq,mtsq,musq,0)
     & +6._dp*Bdiff(s23,zip,mtsq)+12._dp)
      return
      end
