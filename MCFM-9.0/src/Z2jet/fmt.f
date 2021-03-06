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
 
      function fmt(s12,s34,s56)
      implicit none
      include 'types.f'
      complex(dp):: fmt
      include 'masses.f'
      real(dp):: s12,s34,s56
! The full BDK expression keeps terms up to 1/mt^4 (II.15) but here we
! we will only use up to 1/mt^2, as done elsewhere
!      fmt=((1._dp+(2._dp*s34+s12+s56)/15._dp/mt**2)/(24._dp*mt**2))
      fmt=((1._dp)/(24._dp*mt**2))
      return
      end

