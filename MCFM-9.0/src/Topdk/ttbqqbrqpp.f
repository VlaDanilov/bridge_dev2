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
 
      function ttbqqbrqpp(k1,k2,k3,k4,k5,k6,k7)
      implicit none
      include 'types.f'
      complex(dp):: ttbqqbrqpp
      
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'masses.f'
      integer:: k1,k2,k3,k4,k5,k6,k7
      real(dp):: s3459,s6789,mtsq
      s3459=s(k4,k3)+s(k5,k3)
      s6789=s(k3,k6)+s(k3,k7)
      mtsq=mt**2
      ttbqqbrqpp =  + xn**(-1)*s6789**(-1) * ( 1/(za(k1,k2))/(za(k5,k3)
     &    )/(zb(k1,k2))*za(k1,k5)*za(k5,k6)*za(k6,k7)*zb(k2,k6)*zb(k4,
     &    k5)*zb(k6,k3) + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,
     &    k5)*za(k5,k7)*za(k6,k7)*zb(k2,k7)*zb(k4,k5)*zb(k6,k3) + 1/(
     &    za(k1,k2))/(zb(k1,k2))*za(k1,k5)*za(k6,k7)*zb(k2,k3)*zb(k4,k5
     &    )*zb(k6,k3) )
      ttbqqbrqpp = ttbqqbrqpp + xn**(-1)*mtsq*s6789**(-1) * ( 1/(za(k1,
     &    k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k5,k7)*zb(k2,k3)*
     &    zb(k4,k5) + 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*
     &    za(k6,k7)*zb(k2,k4)*zb(k6,k3) - 1/(za(k1,k2))/(za(k5,k3))/(
     &    zb(k1,k2))*za(k1,k6)*za(k5,k7)*zb(k2,k4)*zb(k6,k3) - 1/(za(k1
     &    ,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k7)*za(k5,k7)*zb(k2,k4)*
     &    zb(k7,k3) )
      ttbqqbrqpp = ttbqqbrqpp + xn**(-1)*mtsq*s3459**(-1) * (  - 1/(za(
     &    k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k5)*za(k6,k7)*zb(k2,k6)
     &    *zb(k4,k3) - 1/(za(k1,k2))/(za(k5,k3))/(zb(k1,k2))*za(k1,k7)*
     &    za(k4,k5)*zb(k2,k4)*zb(k4,k3) + 1/(za(k1,k2))/(zb(k1,k2))*za(
     &    k1,k7)*zb(k2,k3)*zb(k4,k3) )

      return
      end
