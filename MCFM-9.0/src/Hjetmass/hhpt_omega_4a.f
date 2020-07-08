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
      w(2)=Log(1 - u4)
      w(3)=Log(mt2*u**(-1))
      w(4)=Log(u4)
      w(5)=1/( - 2*s + 6*s*u4 - 6*s*u4**2 + 2*s*u4**3)
      w(6)=im*Pi
      w(7)=2*w(6)
      w(8)=w(7)*w(2)
      w(9)=w(2)**2
      w(10)=Pi**2
      w(11)=2*w(10) + w(8) + w(9)
      w(11)=w(1)*w(11)
      w(12)=Pi*w(1)
      w(13)=im*w(1)
      w(14)=w(2)*w(13)
      w(12)=w(12) + w(14)
      w(14)=2*im
      w(12)=w(12)*w(14)
      w(15)=w(4)*w(1)
      w(12)=w(12) + w(15)
      w(12)=w(4)*w(12)
      w(16)=w(2) - w(6)
      w(16)=w(1)*w(16)
      w(16)=w(15) + w(16)
      w(17)=w(3)*w(1)
      w(16)=2*w(16) - 3*w(17)
      w(16)=w(3)*w(16)
      w(11)=w(16) + w(12) + w(11)
      w(12)=2*mt2
      w(11)=w(11)*w(12)
      w(16)=4*w(6)
      w(17)=w(2) - 2
      w(18)= - w(16) - w(17)
      w(19)=7*w(3)
      w(20)=5*w(4)
      w(18)= - w(19) + 3*w(18) + w(20)
      w(18)=w(3)*w(18)
      w(17)=3*w(17)
      w(21)= - w(6)*w(17)
      w(22)=w(9) + 12
      w(23)=5*w(6)
      w(24)=w(23) - w(4)
      w(25)=w(24) + 2*w(2)
      w(25)=w(4)*w(25)
      w(18)=w(18) + w(25) + w(21) + 5*w(10) - w(22)
      w(21)=4*mt2
      w(18)=w(18)*w(21)
      w(25)=w(2) - 7
      w(26)=w(25)*w(6)
      w(27)=w(25) + 6*w(3)
      w(28)= - w(4) + 12*w(6) + w(27)
      w(28)=w(3)*w(28)
      w(29)=w(2) + 12
      w(30)= - 1 - w(6)
      w(30)=w(4)*w(30)
      w(26)=w(28) + w(30) + w(26) - 6*w(10) + w(29)
      w(26)=w(26)*w(12)
      w(28)=w(3) - 1
      w(30)= - w(7) - w(28)
      w(30)=w(3)*w(30)
      w(30)=w(30) + w(6) - 2 + w(10)
      w(30)=u4*w(30)*w(21)
      w(16)=w(16) + 12 - w(10)
      w(16)=s*w(16)
      w(31)=2*s
      w(32)=2 + w(6)
      w(32)=w(32)*w(31)
      w(33)=w(3)*s
      w(32)=w(32) + w(33)
      w(32)=w(3)*w(32)
      w(16)=w(30) + w(26) + w(16) + w(32)
      w(16)=u4*w(16)
      w(26)= - w(29)*w(7)
      w(30)=w(9) + 72
      w(32)=w(2) + w(6)
      w(32)=2*w(32) - w(4)
      w(32)=w(4)*w(32)
      w(26)=w(32) + w(26) + 4*w(10) - w(30)
      w(26)=s*w(26)
      w(24)= - w(29) - w(24)
      w(24)=w(24)*w(31)
      w(32)=5*w(33)
      w(24)=w(24) - w(32)
      w(24)=w(3)*w(24)
      w(16)=2*w(16) + w(18) + w(26) + w(24)
      w(16)=u4*w(16)
      w(17)=w(17) + w(23)
      w(18)=12*w(4)
      w(19)=w(19) - w(18) + w(17)
      w(19)=w(3)*w(19)
      w(17)=w(20) - w(17)
      w(17)=w(4)*w(17)
      w(17)=w(19) + w(17) + w(8) - w(10) + w(22)
      w(12)=w(17)*w(12)
      w(17)=w(2) + 4
      w(19)=2*w(4)
      w(22)=w(19) + w(6) - w(17)
      w(22)=w(22)*w(19)
      w(17)=w(17)*w(7)
      w(9)=w(22) + w(17) - w(10) + 48 + w(9)
      w(9)=s*w(9)
      w(7)=w(7) + 8 - w(19) + w(2)
      w(7)=w(7)*w(31)
      w(7)=w(7) + 3*w(33)
      w(7)=w(3)*w(7)
      w(7)=w(12) + w(9) + w(7)
      w(7)=2*w(7) + w(16)
      w(7)=u4*w(7)
      w(9)=w(18) - w(6) - w(27)
      w(9)=w(3)*w(9)
      w(10)=w(29) + w(6)
      w(6)= - 6*w(4) + w(25) + w(6)
      w(6)=w(4)*w(6)
      w(6)=w(9) + w(6) - w(10)
      w(6)=w(6)*w(21)
      w(9)=2*w(10) - w(20)
      w(9)=w(4)*w(9)
      w(8)=w(9) - w(8) - w(30)
      w(8)=s*w(8)
      w(9)=w(20) - w(10)
      w(9)=w(9)*w(31)
      w(9)=w(9) - w(32)
      w(9)=w(3)*w(9)
      w(6)=w(7) + w(6) + w(8) + w(9)
      w(6)=u4*w(6)
      w(7)=2 - w(4)
      w(7)=w(7)*w(31)
      w(7)=w(7) + w(33)
      w(7)=w(3)*w(7)
      w(8)= - w(19) + w(28)
      w(8)=w(3)*w(8)
      w(9)=w(4) + 1
      w(9)=w(9)*w(4)
      w(8)=w(8) + 2 + w(9)
      w(8)=w(8)*w(21)
      w(10)= - 4 + w(4)
      w(10)=w(4)*w(10)
      w(10)=12 + w(10)
      w(10)=s*w(10)
      w(7)=w(8) + w(10) + w(7)
      w(6)=2*w(7) + w(6)
      w(6)=w(5)*w(6)
      w(7)=w(13)*w(28)
      w(8)= - w(14)*w(15)
      w(7)=w(8) + w(7)
      w(7)=w(3)*w(7)
      w(8)=w(13)*w(9)
      w(9)=w(1)*w(14)
      w(7)=w(7) + w(9) + w(8)
      w(7)=w(7)*w(21)
      w(8)=w(4)*im
      w(9)=w(14) - w(8)
      w(10)=w(3)*im
      w(9)=2*w(9) + w(10)
      w(9)=w(3)*w(9)
      w(8)= - 4*im + w(8)
      w(8)=w(4)*w(8)
      w(7)=w(7) + w(9) + 12*im + w(8)


      omega_ggg_ppp_1l_4a = w(11)

      omega_ggg_pmp_1l_4a = w(6)

      omega_qag_mpp_1l_4a = w(7)