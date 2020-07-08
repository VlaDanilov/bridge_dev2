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
 
      subroutine finitemtcorr(rescaling)
      implicit none
      include 'types.f'
      
      real(dp):: rescaling,tn
      complex(dp):: ftn
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      
      tn = 4._dp*(mt/hmass)**2
      if (tn < (1._dp)) then
         ftn = 0.5_dp*(log((1._dp+sqrt(1._dp-tn))
     &        /(1._dp-sqrt(1._dp-tn)))-im*pi)**2
      else
         ftn = -2._dp*(asin(1.0_dp/sqrt(tn)))**2
      endif
      rescaling=abs((3._dp*tn/4._dp)*(2._dp+(tn-1._dp)*ftn))**2
      
      return
      end
      
