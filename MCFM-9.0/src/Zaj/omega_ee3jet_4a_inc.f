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
      w(2)=OneMz**(-1)
      w(3)=y**(-1)
      w(4)=log( - 1 + xv**(-1))
      w(5)=log(1 + xu*xv**(-1))
      w(6)=log(1 - xu)
      w(7)=log(1 - xv - xu)
      w(8)=log(1 + 1/( - 1 + xv)*xu)
      w(9)=log(xu)
      w(10)=log(xv)
      w(11)=li2(1 - xv**(-1))
      w(12)=li2(1 - xv**(-1) + xu*xv**(-1))
      w(13)=li2( - xu*xv**(-1))
      w(14)=yPz**(-1)
      w(15)=w(10) - w(7)
      w(16)=pi*im
      w(17)=w(15) - w(16)
      w(17)=w(17)*w(6)
      w(17)=w(17) - w(12)
      w(18)=w(4) + w(8)
      w(19)=w(18) + w(16)
      w(20)=w(19)*w(9)
      w(20)=w(20) + w(11)
      w(15)=w(18) + w(15)
      w(21)=w(15) - 3._dp/4._dp
      w(21)=w(21)*w(10)
      w(18)=1._dp/2._dp*w(18)
      w(22)=w(16) + 3._dp/2._dp
      w(22)=w(18)*w(22)
      w(20)= - w(21) + w(22) + 1._dp/2._dp*w(20) + w(17)
      w(21)=im**2
      w(21)=w(21) + 1._dp/6._dp
      w(21)=w(21)*pi
      w(22)=3._dp/2._dp*im
      w(21)=w(21) + w(22)
      w(21)=w(21)*pi
      w(23)=w(2)*w(10)
      w(24)=1._dp/2._dp*w(23)
      w(25)=w(10) + w(24) - 1._dp/2._dp
      w(25)=w(25)*w(2)
      w(26)=w(13) - 7._dp/2._dp
      w(27)=w(25) - w(21) + w(26)
      w(28)=1._dp/2._dp*w(3)
      w(29)=w(23) - w(10)
      w(30)= - w(29)*w(28)
      w(27)=w(30) + w(20) - 1._dp/2._dp*w(27)
      w(27)=w(1)*w(27)
      w(22)= - w(22) + 1._dp/3._dp*pi
      w(30)=1._dp/2._dp*pi
      w(22)=w(22)*w(30)
      w(30)=w(16) + w(9)
      w(31)=w(10) - w(30)
      w(32)=w(31)*w(10)
      w(33)=Pi**2
      w(31)=w(31)*w(5)
      w(22)= - w(22) + 3._dp/4._dp*w(9) + w(32) + w(31) + w(11) + 1._dp/6.D
     & 0*w(33)
      w(25)=w(25) + w(26) - w(22)
      w(26)=1._dp/2._dp*N
      w(25)=w(25)*w(26)
      w(32)=w(15) - 1._dp/2._dp
      w(32)=w(32)*w(10)
      w(16)=w(18) + 1._dp/2._dp*w(16)
      w(17)=w(11) + w(17)
      w(18)=w(32) - w(16) - w(17)
      w(32)=w(1)*w(3)
      w(33)= - w(18)*w(32)
      w(15)=w(15)*w(10)
      w(15)=w(17) - w(15)
      w(17)=w(15)*w(3)**2
      w(34)=1._dp/2._dp*w(1)
      w(35)=z*w(34)*w(17)
      w(33)=w(33) + w(35)
      w(33)=z*w(33)
      w(36)=Pi*im
      w(30)= - w(36) - w(10) + 1._dp/2._dp*w(30)
      w(36)=1._dp/2._dp*beta0
      w(30)=w(30)*w(36)
      w(36)=w(34)*w(31)
      w(25)=w(33) + w(25) + w(36) - w(30) + w(27)
      w(16)= - w(16) - w(15)
      w(16)=w(3)*w(16)
      w(27)=w(13) - 3
      w(21)=w(21) - w(27)
      w(16)=w(16) + 1._dp/4._dp*w(23) + 1._dp/2._dp*w(21) + w(20)
      w(16)=w(1)*w(16)
      w(20)= - w(24) + w(27) - w(22)
      w(20)=w(20)*w(26)
      w(15)= - w(15)*w(28)
      w(15)=w(15) - w(18)
      w(15)=w(15)*w(32)
      w(18)=w(19)*w(14)
      w(21)=w(1)*w(18)
      w(22)=w(21) + w(1)
      w(22)=w(22)*w(14)
      w(15)=w(35) + w(15) + 1._dp/2._dp*w(22)
      w(15)=z*w(15)
      w(18)= - w(18) + w(31)
      w(18)=w(34)*w(18)
      w(15)=w(15) + w(20) - w(30) + w(16) + w(18)
      w(16)= - w(1)*w(17)
      w(16)=w(16) - w(22)
      w(16)=z*w(16)
      w(17)= - w(19) + w(29)
      w(17)=w(3)*w(17)
      w(18)=w(29) - 1
      w(18)=w(18)*w(2)
      w(18)=w(18) + 1
      w(17)=1._dp/2._dp*w(18) + w(17)
      w(17)=w(1)*w(17)
      w(18)= - w(18)*w(26)
      w(16)=w(16) + w(21) + w(17) + w(18)
      w(16)=1._dp/2._dp*w(16)


       c_alphai(1) = w(25)

       c_betai(1) = w(15)

       c_gammai(1) = w(16)

       c_alphai(0) =  1._dp

       c_betai(0) =  1._dp

       c_gammai(0) =  0
