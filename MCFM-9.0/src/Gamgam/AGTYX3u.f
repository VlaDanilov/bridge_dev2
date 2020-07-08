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
 
      function AGTYX3u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. B.14
      include 'types.f'
      real(dp):: AGTYX3u
      include 'constants.f'
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYX3u= (1._dp/ 18*Lx**2+ (-2._dp/ 9*Ly+2._dp/ 9*Lu )*Lx
     & -4._dp/9*Ly*Lu+1._dp/18*pisq+2._dp/9*Lu**2+2._dp/9*Ly**2)
     & *(t**2+s**2)/(s*t)
      return
      end
