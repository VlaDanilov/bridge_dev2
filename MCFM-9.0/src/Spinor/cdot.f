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
 
      function cdot(e1,e2)
      implicit none
      include 'types.f'
      complex(dp):: cdot
      
      complex(dp):: e1(4),e2(4)
      cdot=e1(4)*e2(4)-e1(1)*e2(1)-e1(2)*e2(2)-e1(3)*e2(3)
      return
      end
