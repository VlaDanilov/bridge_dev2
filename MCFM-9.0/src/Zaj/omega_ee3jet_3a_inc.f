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
 
      !complex(dp) :: w(36)

      w(1)=N**(-1)
      w(2)=y**(-1)
      w(3)=OneMz**(-1)
      w(4)=log( - 1 + 1/(xv + xu))
      w(5)=log(1 + xu*xv**(-1))
      w(6)=log(1 - xu)
      w(7)=log(1 - xv - xu)
      w(8)=log(xu)
      w(9)=log(xv)
      w(10)=li2(1 - xv**(-1))
      w(11)=li2(1 - xv**(-1) + xu*xv**(-1))
      w(12)=li2( - xu*xv**(-1))
      w(13)=yPz**(-1)
      w(14)=im*pi
      w(15)=w(9) - w(7)
      w(16)=w(15) - w(14)
      w(16)=w(16)*w(6)
      w(16)=w(16) - w(11)
      w(17)=w(8) + 3._dp/2._dp
      w(18)= - w(9) + 1._dp/2._dp*w(14)
      w(17)=w(18) + 1._dp/2._dp*w(17)
      w(17)=w(17)*w(4)
      w(19)=w(9) - 3._dp/2._dp
      w(20)=1._dp/2._dp*w(5)
      w(19)=w(19)*w(20)
      w(21)=pi**2
      w(22)=w(21)*im
      w(23)=pi*w(8)
      w(23)=w(22) + w(23)
      w(24)=1._dp/2._dp*im
      w(24)=w(23)*w(24)
      w(25)=w(15) - 3._dp/4._dp
      w(25)=w(25)*w(9)
      w(17)=w(17) - w(19) - w(25) + w(24) + 1._dp/2._dp*w(10) + w(16)
      w(19)=w(23)*im
      w(21)=1._dp/6._dp*w(21)
      w(16)=w(21) - w(12) + w(16)
      w(19)=w(19) + w(16)
      w(23)=w(9) - w(8)
      w(24)=w(23) - w(14)
      w(25)=w(24) - 1._dp/2._dp
      w(26)=w(25)*w(4)
      w(20)=w(26) - w(20)
      w(26)=w(15) - 1._dp/2._dp
      w(26)=w(26)*w(9)
      w(27)=1._dp/2._dp*w(8)
      w(26)=w(26) + w(27) - w(19) + w(20)
      w(28)= - z*w(26)
      w(29)=w(24)*w(3)
      w(30)=w(29) - w(24)
      w(31)=w(24)*w(4)
      w(15)=w(15)*w(9)
      w(19)=w(19) - w(31) - w(15)
      w(31)=1._dp/2._dp*w(2)
      w(31)=w(31)*w(19)
      w(32)=z**2*w(31)
      w(28)=w(32) - 1._dp/2._dp*w(30) + w(28)
      w(28)=w(2)*w(28)
      w(32)=1._dp/2._dp*w(3)
      w(33)=w(32)*w(24)
      w(25)=w(33) + w(25)
      w(32)= - w(25)*w(32)
      w(34)= - w(21) + 3._dp/2._dp*w(8)
      w(35)=w(12) - 7._dp/2._dp
      w(36)= - w(35) - w(34)
      w(28)=w(28) + w(32) + 1._dp/2._dp*w(36) + w(17)
      w(28)=w(1)*w(28)
      w(32)=w(9) + 3._dp/4._dp
      w(32)=w(32)*w(14)
      w(24)=w(24)*w(5)
      w(36)=w(23)*w(9)
      w(21)=w(24) - w(32) + w(10) - 3._dp/4._dp*w(8) + w(36) - w(21)
      w(24)=w(3)*w(25)
      w(24)=w(24) + w(35) - w(21)
      w(24)=N*w(24)
      w(25)=beta0*im
      w(32)=Pi*N
      w(25)= - w(25) + 1._dp/6._dp*w(32)
      w(25)=w(25)*Pi
      w(18)=w(18) + w(27)
      w(18)=w(18)*beta0
      w(18)=w(25) + w(18)
      w(24)=w(24) - w(18)
      w(24)=1._dp/2._dp*w(24) + w(28)
      w(25)=z - 1
      w(25)=w(25)*w(31)
      w(25)=w(25) - w(26)
      w(25)=z*w(25)
      w(26)=w(8) + 1._dp/2._dp
      w(26)= - pi*w(26)
      w(22)=w(26) - w(22)
      w(22)=im*w(22)
      w(15)=w(22) + w(15) + w(25) - w(16) + w(20)
      w(15)=w(2)*w(15)
      w(16)=w(5) + w(4)
      w(14)=w(16) + w(14)
      w(20)=w(14)*w(13)
      w(20)=w(20) + 1
      w(20)=w(20)*z
      w(14)=w(20) - w(14)
      w(14)=w(14)*w(13)
      w(20)=w(12) - 3
      w(22)=w(14) - w(20) - w(34)
      w(15)=w(15) + 1._dp/4._dp*w(29) + w(17) + 1._dp/2._dp*w(22)
      w(15)=w(1)*w(15)
      w(17)= - w(33) + w(20) - w(21)
      w(17)=N*w(17)
      w(17)=w(17) - w(18)
      w(15)=1._dp/2._dp*w(17) + w(15)
      w(17)= - w(2)*z*w(19)
      w(16)=w(17) + w(29) - w(23) - w(16)
      w(16)=w(2)*w(16)
      w(17)=w(30) - 1
      w(17)=w(17)*w(3)
      w(17)=w(17) + 1
      w(17)=1._dp/2._dp*w(17)
      w(14)=w(16) + w(17) - w(14)
      w(14)=w(1)*w(14)
      w(16)= - N*w(17)
      w(14)=w(16) + w(14)
      w(14)=1._dp/2._dp*w(14)


       c_alphai(1) = w(24)

       c_betai(1) = w(15)

       c_gammai(1) = w(14)

       c_alphai(0) =  1._dp

       c_betai(0) =  1._dp

       c_gammai(0) =  0
