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
 
      function Gamma0(mt,besq,omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0
C--   Author: John M. Campbell and R.K. Ellis, January 2012  
C--   Taken from formula (2) of
C--   Fermilab-PUB-12-078-T
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      real(dp):: mt,omsq,besq,Gammainfty,f,P3b
      Gammainfty=GF*mt**3/(8._dp*rt2*pi)
      P3b=0.5_dp*sqrt(1._dp+omsq**2+besq**2-2._dp*(omsq+besq+omsq*besq))
      f=(1._dp-besq)**2+omsq*(1._dp+besq)-2._dp*omsq**2
      Gamma0=Gammainfty*2._dp*P3b*f
      return
      end

