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
 
      !complex(dp) w(32)

      w(1)=s**(-1)
      w(2)=Log(mt2*t**(-1))
      w(3)=Hr1(0)
      w(4)=Hr1(1)
      w(5)=Hr2(0,0)
      w(6)=Hr2(0,1)
      w(7)=Hr2(1,0)
      w(8)=Hr2(1,1)
      w(9)=1/( - 2 + 2*u3)
      w(10)=1/( - 1 + u3)
      w(11)=1/( - s + s*u3)
      w(12)=1/(s - 2*s*u3 + s*u3**2)
      w(13)=im*Pi
      w(14)=w(4) + w(3)
      w(15)= - w(2) + w(14)
      w(15)=w(15)*w(13)
      w(16)=w(7) + w(6)
      w(17)=Pi**2
      w(15)= - w(16) + w(17) + w(15)
      w(18)=2*w(2)
      w(14)=w(18)*w(14)
      w(19)=2*w(8)
      w(20)=w(2)**2
      w(21)=2*w(5)
      w(14)=w(21) + w(19) - 3*w(20) + w(14) + 2*w(15)
      w(14)=2*w(1)*w(14)
      w(15)=u3 - 2
      w(22)=w(15)*u3
      w(23)=w(22) + 1
      w(24)=u3**2
      w(17)=w(23)*w(17)*w(24)
      w(23)=w(15)*w(24)
      w(25)=w(2)*u3
      w(26)= - 1 - u3
      w(26)=w(26)*w(25)
      w(23)=4*w(23) + w(26)
      w(23)=w(4)*w(23)
      w(26)=3 + w(2)
      w(27)=u3 - 1
      w(28)=w(27)*u3
      w(29)=w(28) + 1
      w(26)=w(29)*w(26)
      w(29)=4*u3
      w(30)=w(15)*w(29)
      w(30)=5 + w(30)
      w(30)=u3*w(30)
      w(30)=1 + w(30)
      w(30)=w(8)*u3*w(30)
      w(31)=w(25)*w(27)
      w(32)= - w(4)*w(28)
      w(32)=w(31) + w(32)
      w(32)=w(32)*w(13)
      w(23)=w(32) + 2*w(17) + w(23) + 4*w(26) + w(30)
      w(23)=w(10)*w(23)
      w(26)=w(24) - 1
      w(30)=u3*w(13)
      w(15)= - w(2)*w(15)
      w(15)= - w(30) - 4*w(26) + w(15)
      w(15)=w(3)*w(15)
      w(27)=w(29)*w(27)
      w(27)=w(27) + 1
      w(27)=w(27)*u3
      w(29)= - w(16)*w(27)
      w(32)=2 + 3*w(28)
      w(20)=w(9)*w(32)*w(20)
      w(27)= - 2 + w(27)
      w(27)=w(5)*w(27)
      w(15)=w(23) + w(27) + w(20) + w(15) + w(29)
      w(20)=w(22) - 1
      w(19)=w(20)*w(24)*w(19)
      w(23)=u3 - 3
      w(23)=w(23)*u3
      w(27)= - 2 - w(23)
      w(27)=w(4)*w(27)*w(24)
      w(29)= - 4 - 3*w(22)
      w(29)=u3*w(29)
      w(29)=1 + w(29)
      w(29)=u3*w(29)
      w(27)=w(27) + w(29) - w(31)
      w(13)=w(27)*w(13)
      w(27)= - 2 - 5*w(22)
      w(27)=u3*w(27)
      w(27)= - 3 + w(27)
      w(27)=u3*w(27)
      w(22)= - 3 + w(22)
      w(22)=u3*w(22)
      w(22)=4 + w(22)
      w(22)=u3*w(22)
      w(22)= - 2 + w(22)
      w(22)=w(2)*w(22)
      w(22)=w(22) + 2 + w(27)
      w(22)=w(2)*w(22)
      w(23)=w(23) - 1
      w(23)=w(23)*u3
      w(27)=1 - w(23)
      w(25)=w(27)*w(25)
      w(27)=w(24) - 4
      w(27)=w(27)*u3
      w(29)=1 + w(27)
      w(29)=u3*w(29)
      w(25)=w(29) + w(25)
      w(25)=w(4)*w(25)
      w(20)=u3*w(20)
      w(20)=2 + w(20)
      w(20)=u3*w(20)
      w(20)= - 1 + w(20)
      w(13)=w(13) + w(17) + w(25) + w(19) + 4*w(20) + w(22)
      w(13)=w(12)*w(13)
      w(17)= - w(26)*w(30)
      w(19)= - 4 - w(27)
      w(19)=w(2)*w(19)
      w(17)=w(17) + w(19) + 2 + w(23)
      w(17)=w(3)*w(17)
      w(19)= - 1 + 2*w(28)
      w(16)= - u3*w(19)*w(16)
      w(16)=w(16) + w(17)
      w(16)=w(11)*w(16)
      w(17)= - 2 + w(24)
      w(17)=w(1)*w(17)*w(21)
      w(13)=w(13) + w(17) + w(16)
      w(13)=2*w(13)
      w(16)=4 + w(2)
      w(16)=w(2)*w(16)
      w(17)= - 2 - w(2)
      w(17)=w(3)*w(17)
      w(16)=w(21) + 2*w(17) + 12 + w(16)
      w(16)=im*w(16)
      w(17)= - 1 + w(2)
      w(17)=w(2)*w(17)
      w(18)=1 - w(18)
      w(18)=w(3)*w(18)
      w(17)=w(21) + w(18) + 2 + w(17)
      w(17)=4*w(1)*im*w(17)


      omega_ggg_ppp_1l_3a_mt0 =  0

      omega_ggg_pmp_1l_3a_mt0 = w(15)

      omega_qag_mpp_1l_3a_mt0 = w(16)

      omega_ggg_ppp_1l_3a_mt2 = w(14)

      omega_ggg_pmp_1l_3a_mt2 = w(13)

      omega_qag_mpp_1l_3a_mt2 = w(17)
