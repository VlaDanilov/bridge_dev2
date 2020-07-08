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
 
      function F1anom(s12,s45,mt2,musq)
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      complex(dp)::F1anom
      complex(dp)::lnrat
      real(dp)::s12,s45,mt2,musq

      if (mt2 .eq. zero) then
      F1anom=one/two/(s45-s12)*(1.d0+s45/(s45-s12)*lnrat(-s12,-s45))
      return
      else
      F1anom=one/two/(s45-s12)
     & *(one+two*mt2*qlI3(s12,0.d0,s45,mt2,mt2,mt2,musq,0)
     & +(s45/(s45-s12))*(qlI2(s45,mt2,mt2,musq,0)-qlI2(s12,mt2,mt2,musq,0)))
      endif
            
      return
      end

