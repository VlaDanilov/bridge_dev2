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
 
      subroutine breitw1(msq,mminsq,mmaxsq,rmass,rwidth,wt)       
      implicit none
      include 'types.f'
      
c---- Given a mass-squared msq generate a weight wt 
c---- breit-wigner should still be included in the matrix element
c      real(dp):: bw
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zerowidth.f'
      real(dp):: mminsq,mmaxsq,rmass,rwidth,msq,wt
      real(dp):: almin,almax,tanal

      if (zerowidth) then
          tanal=0._dp
          almax=+pi/two
          almin=-pi/two
      else
          almin=atan((mminsq-rmass**2)/rmass/rwidth)
          almax=atan((mmaxsq-rmass**2)/rmass/rwidth)
          tanal=(rmass**2-msq)/(rmass*rwidth)
      endif

      wt=(almax-almin)*rmass*rwidth*(1._dp+tanal**2)
      return
      end

