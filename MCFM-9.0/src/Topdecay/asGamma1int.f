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
 
      function asGamma1int(omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1int
      
C--   Author R.K. Ellis April 2012
C--   Integrand for NLO width with W-offshell 
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      real(dp):: omsq,mt,xi,ga,besq,asGamma1
      common/transfer/mt,besq,xi,ga
!$omp threadprivate(/transfer/)
      asGamma1int=ga*xi/pi
     & /((1._dp-xi*omsq)**2+ga**2)*asGamma1(mt,besq,omsq)
      return
      end
