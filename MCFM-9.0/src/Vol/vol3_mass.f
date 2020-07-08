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
 
      function vol3_mass(decaymass,s)
      implicit none
      include 'types.f'
      real(dp):: vol3_mass
      
      real(dp):: decaymass,b0,bp,bm,s,ro,vol
c  Volume of phase space for 2-->p3(decaymass)+p4(decaymass)+p5(0)      
      if (decaymass == 0._dp) then
      vol3_mass=1._dp
      else
      ro=4._dp*decaymass**2/s
      b0=sqrt(1._dp-ro)
      bp=0.5_dp*(1._dp+b0)      
      bm=0.5_dp*(1._dp-b0)      

      vol3_mass=(3._dp-b0**2)/2._dp*b0
     & -(1._dp-b0**2)*(3._dp+b0**2)/4._dp*log(bp/bm)
      endif
      vol3_mass=vol3_mass*vol(s,3)
      return
      end
      
