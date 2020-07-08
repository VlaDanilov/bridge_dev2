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
 
      !complex(dp) w(36)

      w(1)=s**(-1)
      w(2)=Log(1 - u2)
      w(3)=Log(mt2*s**(-1))
      w(4)=Log(u2)
      w(5)=1/( - 2*s*u2**2 + 6*s*u2**3 - 6*s*u2**4 + 2*s*u2**5)
      w(6)=w(2)**2
      w(7)=w(4)**2
      w(7)=w(6) + w(7)
      w(7)=w(1)*w(7)
      w(8)=im**2
      w(9)=w(8)*w(1)
      w(10)=2*w(4)
      w(11)=w(2)*w(10)*w(9)
      w(12)=w(2) + w(4)
      w(12)=w(1)*w(12)
      w(13)=im*w(12)
      w(14)=Pi*w(1)
      w(13)=w(13) + w(14)
      w(14)=2*Pi
      w(13)=w(13)*w(14)
      w(15)=Pi*im
      w(16)=w(15)*w(1)
      w(12)= - w(16) + w(12)
      w(17)=w(3)*w(1)
      w(12)=2*w(12) - 3*w(17)
      w(12)=w(3)*w(12)
      w(7)=w(12) + w(13) + w(11) + w(7)
      w(11)=2*mt2
      w(7)=w(7)*w(11)
      w(12)=w(2) + 12
      w(13)=w(12) - w(4)
      w(17)=im*w(13)
      w(17)=w(17) - w(14)
      w(17)=w(17)*w(14)
      w(18)=2*w(2)
      w(19)=w(18) - w(4)
      w(19)=w(19)*w(4)
      w(20)=w(6) + 72
      w(17)=w(17) - w(19) + w(20)
      w(17)=s*w(17)
      w(21)=w(4) - w(2)
      w(22)=w(21) + 7
      w(23)=6*w(3)
      w(24)=12*w(15)
      w(25)=w(23) + w(24) - w(22)
      w(25)=w(3)*w(25)
      w(22)= - im*w(22)
      w(22)=w(22) - 6*Pi
      w(22)=Pi*w(22)
      w(22)=w(25) + w(22) + w(13)
      w(25)=4*mt2
      w(22)=w(22)*w(25)
      w(26)=2*s
      w(27)= - 2 - w(15)
      w(27)=w(27)*w(26)
      w(28)=w(3)*s
      w(27)=w(27) - w(28)
      w(27)=w(3)*w(27)
      w(29)=w(3) - 1
      w(30)=2*w(15)
      w(31)= - w(30) - w(29)
      w(31)=w(3)*w(31)
      w(32)=im + Pi
      w(32)=Pi*w(32)
      w(31)=w(31) - 2 + w(32)
      w(31)=w(31)*w(25)
      w(32)= - 4*im + Pi
      w(32)=Pi*w(32)
      w(32)= - 12 + w(32)
      w(32)=s*w(32)
      w(27)=w(31) + w(32) + w(27)
      w(27)=u2*w(27)
      w(31)=5*w(15)
      w(13)=w(31) + w(13)
      w(13)=w(13)*w(26)
      w(32)=5*w(28)
      w(13)=w(13) + w(32)
      w(13)=w(3)*w(13)
      w(13)=2*w(27) + w(22) + w(17) + w(13)
      w(13)=u2*w(13)
      w(17)=w(2) - 2
      w(22)=5*w(4)
      w(27)= - w(22) + 3*w(17)
      w(33)=7*w(3)
      w(24)= - w(33) - w(24) - w(27)
      w(24)=w(3)*w(24)
      w(34)= - im*w(27)
      w(34)=w(34) + 5*Pi
      w(34)=Pi*w(34)
      w(35)=w(6) + 12
      w(19)=w(24) + w(34) + w(19) - w(35)
      w(19)=w(19)*w(11)
      w(24)=w(2) + 4
      w(34)= - w(10) + w(24)
      w(34)=w(34)*w(10)
      w(36)=2*im
      w(24)= - w(4) - w(24)
      w(24)=w(24)*w(36)
      w(24)=w(24) + Pi
      w(24)=Pi*w(24)
      w(6)=w(24) + w(34) - 48 - w(6)
      w(6)=s*w(6)
      w(24)= - w(30) - 8 + w(10) - w(2)
      w(24)=w(24)*w(26)
      w(24)=w(24) - 3*w(28)
      w(24)=w(3)*w(24)
      w(6)=w(19) + w(6) + w(24)
      w(6)=2*w(6) + w(13)
      w(6)=u2*w(6)
      w(13)= - 4*w(4) + w(17)
      w(13)=w(33) + 3*w(13) + w(31)
      w(13)=w(3)*w(13)
      w(17)= - w(4)*w(27)
      w(18)=w(18) - w(22)
      w(18)=im*w(18)
      w(18)=w(18) - Pi
      w(18)=Pi*w(18)
      w(13)=w(13) + w(18) + w(17) + w(35)
      w(13)=w(13)*w(25)
      w(17)= - 2*w(12) + w(22)
      w(17)=w(4)*w(17)
      w(18)= - w(21)*w(30)
      w(17)=w(18) + w(17) + w(20)
      w(17)=s*w(17)
      w(18)=w(15) - w(22) + w(12)
      w(18)=w(18)*w(26)
      w(18)=w(18) + w(32)
      w(18)=w(3)*w(18)
      w(6)=w(6) + w(13) + w(17) + w(18)
      w(6)=u2*w(6)
      w(13)=w(2) - 7
      w(17)= - w(23) - w(15) + 12*w(4) - w(13)
      w(17)=w(3)*w(17)
      w(18)= - 1 + w(4)
      w(18)=w(18)*w(15)
      w(13)= - 6*w(4) + w(13)
      w(13)=w(4)*w(13)
      w(12)=w(17) + w(18) + w(13) - w(12)
      w(11)=w(12)*w(11)
      w(12)= - 2 + w(4)
      w(12)=w(12)*w(26)
      w(12)=w(12) - w(28)
      w(12)=w(3)*w(12)
      w(13)=4 - w(4)
      w(13)=w(4)*w(13)
      w(13)= - 12 + w(13)
      w(13)=s*w(13)
      w(11)=w(11) + w(13) + w(12)
      w(6)=2*w(11) + w(6)
      w(6)=u2*w(6)
      w(10)= - w(10) + w(29)
      w(10)=w(3)*w(10)
      w(11)=1 + w(4)
      w(11)=w(4)*w(11)
      w(10)=w(10) + 2 + w(11)
      w(10)=mt2*w(10)
      w(6)=8*w(10) + w(6)
      w(6)=w(5)*w(6)
      w(10)= - w(9) - w(16)
      w(10)=Pi*w(10)
      w(9)=w(14)*w(9)
      w(11)=w(29)*im*w(1)
      w(9)=w(9) + w(11)
      w(9)=w(3)*w(9)
      w(11)=w(1)*w(36)
      w(9)=w(9) + w(11) + w(10)
      w(9)=w(9)*w(25)
      w(10)=4*w(8) - w(15)
      w(10)=Pi*w(10)
      w(8)=Pi*w(8)
      w(8)=w(36) + w(8)
      w(11)=w(3)*im
      w(8)=2*w(8) + w(11)
      w(8)=w(3)*w(8)
      w(8)=w(9) + w(8) + 12*im + w(10)


      omega_ggg_ppp_1l_2a = w(7)

      omega_ggg_pmp_1l_2a = w(6)

      omega_qag_mpp_1l_2a = w(8)
