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
 
      function AGTYF1u(s,t,Lx,Ly,Lu)
      implicit none
!     Results taken from hep-ph/0201274
!  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
!     Eq. A.18
      include 'types.f'
      real(dp):: AGTYF1u
      include 'constants.f'
      include 'zeta.f'
      real(dp)::s,t,Lx,Ly,Lu

      AGTYF1u=(5._dp/36._dp*Lx**2+ (1._dp/3*Lu
     & -10._dp/27-1._dp/3._dp*Ly )*Lx
     & +1._dp/3._dp*Ly**2+ (20._dp/27-2._dp/3._dp*Lu )*Ly
     & -20._dp/27._dp*Lu-13._dp/108._dp*pisq+1._dp/3._dp*Lu**2)
     & *(t**2+s**2)/(s*t) 
      return
      end
