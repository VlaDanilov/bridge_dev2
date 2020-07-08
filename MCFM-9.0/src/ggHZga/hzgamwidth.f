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
 
      function hzgamwidth(mh)
      implicit none
      include 'types.f'
      real(dp):: hzgamwidth
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'ewcouple.f'
      real(dp):: mh,mhsq
      complex(dp):: f0DDHK
      
      mhsq=mh**2
      
      hzgamwidth=esq*Gf**2*wmass**2*xw/256._dp/pi**5
     & *mh**3*(1._dp-zmass**2/mhsq)**3*abs(f0DDHK(mhsq,zmass**2))**2
      
      return
      end
      
