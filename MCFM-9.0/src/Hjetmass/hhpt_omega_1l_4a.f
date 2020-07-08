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
 
      !complex(dp) w(33)

      w(1)=s**(-1)
      w(2)=Log(mt2*u**(-1))
      w(3)=Hr1(0)
      w(4)=Hr1(1)
      w(5)=Hr2(0,0)
      w(6)=Hr2(0,1)
      w(7)=Hr2(1,0)
      w(8)=Hr2(1,1)
      w(9)=1/( - 2 + 2*u4)
      w(10)=1/( - 1 + u4)
      w(11)=1/( - 1 + 3*u4 - 3*u4**2 + u4**3)
      w(12)=1/(1 - 2*u4 + u4**2)
      w(13)=1/(s - 2*s*u4 + s*u4**2)
      w(14)=im*Pi
      w(15)=w(2) - w(3)
      w(16)=w(4) - w(15)
      w(16)=w(16)*w(14)
      w(17)=w(4)*w(2)
      w(18)=Pi**2
      w(16)=w(16) + w(17) + w(18)
      w(19)=w(7) + w(6)
      w(20)=w(19) - w(8)
      w(21)=2*w(5)
      w(22)= - w(21) + 2*w(20)
      w(23)=3*w(2)
      w(24)=2*w(3)
      w(25)=w(24) - w(23)
      w(25)=w(2)*w(25)
      w(16)=w(25) - w(22) + 2*w(16)
      w(16)=2*w(1)*w(16)
      w(25)= - w(2) - w(4)
      w(25)=w(25)*w(14)
      w(26)= - 4 + w(3)
      w(26)=w(2)*w(26)
      w(20)=w(25) - w(17) + w(26) + w(20)
      w(20)=w(10)*w(20)
      w(25)=2*w(2)
      w(26)=w(2)*w(14)
      w(26)=w(25) + w(26)
      w(26)=w(10)*w(26)
      w(27)=4*w(3)
      w(28)=w(14)*w(3)
      w(29)=4*w(5) - w(27) + w(28)
      w(29)=w(11)*w(29)
      w(30)=w(9)*w(2)**2
      w(26)=w(29) + w(26) + w(30)
      w(29)= - 8*w(14) - 24 + w(18)
      w(29)=w(12)*w(29)
      w(31)=4*w(14) + 12 - w(18)
      w(31)=w(12)*w(31)
      w(32)=w(28) - w(5)
      w(32)=w(11)*w(32)
      w(31)=w(31) + w(32)
      w(31)=u4*w(31)
      w(26)=w(31) + w(29) + 2*w(26)
      w(26)=u4*w(26)
      w(28)= - 5*w(5) + 12*w(3) + w(28)
      w(28)=w(11)*w(28)
      w(20)=w(26) + w(28) + 24*w(12) + w(20) - w(30)
      w(20)=u4*w(20)
      w(26)=w(3) - 2
      w(28)= - w(10)*w(2)*w(26)
      w(29)= - w(24) + w(5)
      w(29)=w(11)*w(29)
      w(28)=w(29) - 6*w(12) + w(28) + w(30)
      w(20)=2*w(28) + w(20)
      w(25)=1 - w(25)
      w(25)=w(25)*w(14)
      w(28)=w(2) - 1
      w(29)= - w(2)*w(28)
      w(25)=w(25) + w(18) - 2 + w(29)
      w(25)=u4*w(25)
      w(29)=w(4) - w(3)
      w(30)=8*w(2) - 5 + w(29)
      w(30)=w(30)*w(14)
      w(31)=4*w(2)
      w(32)=w(31) - 5
      w(33)= - w(3) + w(32)
      w(33)=w(2)*w(33)
      w(25)=2*w(25) + w(30) - 4*w(18) + w(33) + 8 + w(29)
      w(25)=u4*w(25)
      w(27)=w(27) + 1
      w(23)= - w(23) + w(27)
      w(23)=w(2)*w(23)
      w(27)= - 2*w(4) - w(31) + w(27)
      w(27)=w(27)*w(14)
      w(18)=w(25) + w(27) + w(18) + w(23) - 4 + w(29) + w(22)
      w(18)=u4*w(18)
      w(15)=1 + w(15)
      w(14)=w(15)*w(14)
      w(15)= - 8*w(3) + w(32)
      w(15)=w(2)*w(15)
      w(14)=w(18) + 8*w(5) + w(14) + w(4) + w(15) + 8 + 5*w(3) - w(19)
      w(14)=u4*w(14)
      w(15)=w(28) - w(24)
      w(15)=w(15)*w(2)
      w(15)=w(15) + w(3) + 2
      w(18)= - w(21) - w(15)
      w(14)=2*w(18) + w(14)
      w(14)=w(13)*w(14)
      w(17)=w(1)*u4*w(17)
      w(14)=w(17) + w(14)
      w(14)=2*w(14)
      w(17)=3 - w(3)
      w(18)= - 2*w(26) + w(2)
      w(18)=w(2)*w(18)
      w(17)=4*w(17) + w(18)
      w(17)=im*w(17)
      w(18)=w(21)*im
      w(17)=w(17) + w(18)
      w(15)=im*w(15)
      w(15)=w(15) + w(18)
      w(15)=4*w(1)*w(15)


      omega_ggg_ppp_1l_4a_mt0 =  0

      omega_ggg_pmp_1l_4a_mt0 = w(20)

      omega_qag_mpp_1l_4a_mt0 = w(17)

      omega_ggg_ppp_1l_4a_mt2 = w(16)

      omega_ggg_pmp_1l_4a_mt2 = w(14)

      omega_qag_mpp_1l_4a_mt2 = w(15)
