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
 
      function plus(pz,fx,L0,p1,fx1,L01,z,jaco)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: plus
      real(dp), intent(in)::pz,p1,fx,fx1,z,L0,L01,jaco
      if(abs(one-z) < 1.e-5_dp) then
      plus=L01*p1*fx1
      else
      plus=(pz*fx/z-p1*fx1)*L0*jaco+L01*p1*fx1
      endif
      return
      end
