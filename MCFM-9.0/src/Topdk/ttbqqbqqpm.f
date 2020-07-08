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
 
      function ttbqqbqqpm(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'types.f'
      complex(dp):: ttbqqbqqpm
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      real(dp):: s129,mtsq
      s129=s(k1,k2)+s(k1,k3)+s(k2,k3)
      mtsq=mt**2
      ttbqqbqqpm =  + xn**(-1)*s129**(-1) * ( 1/(zb(k1,k3))*za(k5,k3)*
     &    za(k6,k7)*zb(k2,k6)*zb(k4,k5) + 1/(zb(k1,k3))/(zb(k5,k3))*za(
     &    k1,k5)*za(k6,k7)*zb(k1,k5)*zb(k2,k6)*zb(k4,k5) - 1/(zb(k2,k3)
     &    )/(zb(k5,k3))*za(k1,k5)*za(k6,k7)*zb(k2,k5)*zb(k2,k6)*zb(k4,
     &    k5) )
      ttbqqbqqpm = ttbqqbqqpm + xn**(-1)*mtsq*s129**(-1) * (  - 1/(zb(
     &    k1,k3))*za(k7,k3)*zb(k2,k4) - 1/(zb(k1,k3))/(zb(k5,k3))*za(k1
     &    ,k7)*zb(k1,k5)*zb(k2,k4) + 1/(zb(k2,k3))/(zb(k5,k3))*za(k1,k7
     &    )*zb(k2,k4)*zb(k2,k5) )

      return
      end
