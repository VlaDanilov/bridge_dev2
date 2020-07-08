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
 
      function H4prenorm()
      implicit none
      include 'types.f'
      real(dp):: H4prenorm
c--- This function returns the appropriate renormalization factor for
c--- the Higgs + 4 parton amplitudes
c--- it includes:    a) strong coupling renormalization
c---                 b) finite renormalization of Hgg effective coupling
c---                 c) finite renormalization of alpha-s in dred scheme
c---
c--- Note that this function returns zero when checking the
c--- (unrenormalized) results in the EGZ paper
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'b0.f'
      include 'epinv.f'
      include 'scheme.f'
      logical:: CheckEGZ
      common/CheckEGZ/CheckEGZ
!$omp threadprivate(/CheckEGZ/)

      if (CheckEGZ) then
        H4prenorm=0._dp
      else
        H4prenorm=(-4._dp*b0/xn*epinv+11._dp/xn)
        if (scheme == 'dred') then
        H4prenorm=H4prenorm+2._dp/3._dp
        endif
      endif  
      
      return
      end
