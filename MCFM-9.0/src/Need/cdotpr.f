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
 
      function cdotpr(J1,J2)
      implicit none
      include 'types.f'
      complex(dp):: cdotpr
      
      complex(dp):: J1(4),J2(4)
      cdotpr=J1(4)*J2(4)-J1(1)*J2(1)-J1(2)*J2(2)-J1(3)*J2(3)
      return
      end
