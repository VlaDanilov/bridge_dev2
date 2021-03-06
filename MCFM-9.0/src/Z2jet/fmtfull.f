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
 
      function fmtfull(s12,s34,s56) 
        use mod_qcdloop_c
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'scale.f'
      include 'masses.f'
      complex(dp):: fmtfull
      real(dp):: mtsq,s12,s34,s56,Del12,Del34,Del56,Del3
      complex(dp):: bdiff
      include 'cplx.h'

C     Implementation of BDKW Nucl Phys. 513, 3 (1998)
C     Eqn. (B.15) upgraded so that is gives the full mass
C     dependence.
      Del12=s12-s34-s56
      Del34=s34-s12-s56
      Del56=s56-s12-s34
      Del3=s12*Del12+s34*Del34+s56*Del56
      mtsq=mt**2
       
      fmtfull=
     & ((three*s12*s56*s34*Del34/Del3**2-(s12*s56-mtsq*del34)/Del3))
     & *(-qlI3(s12,s34,s56,mtsq,mtsq,mtsq,musq,0))
     & +((three*s56*Del56/Del3**2-half/Del3)*s12)
     & *bdiff(s34,s12,mtsq)
     & +((three*s12*Del12/Del3**2-half/Del3)*s56)
     & *bdiff(s34,s56,mtsq)
     & -cplx2(half*Del34/Del3,zip)
      end

